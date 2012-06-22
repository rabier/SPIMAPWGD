/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gene tree search functions

=============================================================================*/


#include "common.h"
#include "distmatrix.h"
#include "logging.h"
#include "Matrix.h"
#include "model.h"
#include "nj.h"
#include "parsimony.h"
#include "phylogeny.h"
#include "search.h"
#include "top_change.h"
#include "top_prior.h"
#include "treevis.h"

namespace spidir {

//=============================================================================


TreeSet::~TreeSet()
{
    clear();
}

void TreeSet::clear()
{
    for (Set::iterator it=trees.begin();
         it != trees.end(); it++)
    {
        delete [] (int*) *it;
    }

    trees.clear();
}

bool TreeSet::insert(Tree *tree)
{
    int *key2 = new int [tree->nnodes+1];
    tree->hashkey(key2);
    key2[tree->nnodes] = -2; // cap key

    Set::iterator it = trees.find(key2);
    if (it == trees.end()) {
        trees.insert(key2);
        return true;
    } else {
        delete [] key2;
        return false;
    }
}

bool TreeSet::has(Tree *tree)
{
    key.ensureSize(tree->nnodes+1);
    tree->hashkey(key);
    key[tree->nnodes] = -2; // cap key
    return trees.find(key) != trees.end();
}



//=============================================================================

// NNI Proposer

NniProposer::NniProposer(int niter) :
    niter(niter),
    iter(0),
    nodea(NULL),
    nodeb(NULL)
{}

void NniProposer::propose(Tree *tree)
{
    // increase iteration
    iter++;
    // propose new tree
    proposeRandomNni(tree, &nodea, &nodeb);
    performNni(tree, nodea, nodeb);
}

void NniProposer::revert(Tree *tree)
{
  //writeNewickTree(stdout, tree,true);
    // undo topology change
  performNni(tree, nodea, nodeb); 
  //writeNewickTree(stdout, tree,true);

}


//=============================================================================
// SPR Proposer

SprProposer::SprProposer(int niter) :
    NniProposer(niter)
{
}



void SprProposer::propose(Tree *tree)
{   
    // increase iteration
    iter++;
 
    // choose a SPR move
    proposeRandomSpr(tree, &nodea, &nodeb);
    
    // remember sibling of nodea
    const Node *p = nodea->parent;
    nodec = (p->children[0] == nodea) ? p->children[1] : p->children[0];
    
    // perform SPR move
    performSpr(tree, nodea, nodeb);
}

void SprProposer::revert(Tree *tree)
{
    performSpr(tree, nodea, nodec);
}


//=============================================================================
// Mixture of Proposers

void MixProposer::addProposer(TopologyProposer *proposer, float weight)
{
    totalWeight += weight;
    methods.push_back(Method(proposer, weight));
}

void MixProposer::propose(Tree *tree)
{

    // increase iteration
    iter++;

    // randomly choose method
    float choice = frand() * totalWeight;
    float sum = methods[0].second;
    unsigned int i = 0;
    while (i < methods.size()-1 && sum < choice) {
        i++;
        sum += methods[i].second;
    }

    // make proposal
    lastPropose = i;
    methods[i].first->propose(tree);
}

void MixProposer::revert(Tree *tree)
{
    methods[lastPropose].first->revert(tree);
}



//=============================================================================
// SPR Neighborhood Proposer

SprNbrProposer::SprNbrProposer(int niter, int radius) :
    NniProposer(niter),
    radius(radius),
    basetree(NULL),
    reverted(false)
{
}

void SprNbrProposer::propose(Tree *tree)
{
    // ensure the same tree is used for each proposal
    if (!basetree)
        basetree = tree;
    else
        assert(basetree == tree);

    // start a new subtree
    if (iter == 0 || queue.size() == 0 || !reverted) {
        iter = 0;
        pickNewSubtree();
    }

    iter++; // increase iteration
    //printf("\nthis is our non random queue\n");
    list<Node*>::iterator indice=queue.begin();
    for (int i=0; i<queue.size(); i++){
      //printf("\nname of the node %d\n",(*indice)->name);
      indice++;
    }

    //transform queue which is not random by construction
    //to a randomqueue
    Node *node1;
    float randomindex=0;
    list<Node*>::iterator theindex=randomqueue.begin();
    sizequeue=queue.size();    

    for (int i=0; i<sizequeue; i++){
      node1 = queue.front();
      if (randomqueue.size()==0){
	 randomqueue.push_back(node1);
      }else if(randomqueue.size()==1){

	if (frand()>0.5){
	  randomqueue.push_back(node1);}
	else{ 
	  randomqueue.push_front(node1);}
      }
      else{
	randomindex=floor(randomqueue.size()*frand());      
	theindex=randomqueue.begin(); 
	for (int j=0; j<int(randomindex); j++){
	  theindex++;
	}      
	randomqueue.insert(theindex,node1);
      }
       queue.pop_front();
    }
    //we have now our randomqueue
    
    list<Node*>::iterator indices=randomqueue.begin();
    for (int i=0; i<randomqueue.size(); i++){
      //printf("\nname of the node %d\n",(*indices)->name);
      indices++;
    }

    while (randomqueue.size() > 0) {
 
      nodea = randomqueue.front();
      randomqueue.pop_front();      
      // remember sibling of subtree (nodeb)
      const Node *p = subtree->parent;
      nodeb = (p->children[0] == subtree) ? p->children[1] : p->children[0];    
      // perform only valid SPR moves
      // NOTE: the tree may have changed, thus we need to double check
      // whether the Spr is valid.
      if (validSpr(tree, subtree, nodea)) {
	performSpr(tree, subtree, nodea);
	break;
      }
    }

    randomqueue.clear();
    revertsizequeue(tree);//count the number of nodes suitables for a new spr
    assert(tree->assertTree());
}

void SprNbrProposer::revert(Tree *tree)
{

  performSpr(tree, subtree, nodeb);
  assert(tree->assertTree());

}


void SprNbrProposer::revertsizequeue(Tree *tree)
{
    
  Node *a = subtree;

  // find sibling (b) of a
  Node *c = a->parent;
  const int bi = (c->children[0] == a) ? 1 : 0;
  Node *b = c->children[bi];
    
  // uninitialize path distances
  pathdists.clear();
  for (int i=0; i<tree->nnodes; i++)
    pathdists.push_back(-1);
    
  // setup path distances and queue
  pathdists[a->name] = 0;
  pathdists[c->name] = 0;
  pathdists[b->name] = 0;
  queue.clear();
  list<Node*> tmpqueue;
  tmpqueue.push_back(c);
  tmpqueue.push_back(b);

  // traverse tree via depth first traversal
  while (tmpqueue.size() > 0) {
    Node *n = tmpqueue.front();
    tmpqueue.pop_front();
        
    // do not traverse beyond radius
    if (pathdists[n->name] >= radius)
      continue;

    // queue only valid new branch points:
    // n must not be root, a, descendant of a, c (parent of a), or  
    // b (sibling of a)
    if (n->parent && n != subtree && n != b && n != c) {
      queue.push_back(n);
    }

    // queue up unvisited neighboring edges
    Node *w = n->parent;

    if (w && pathdists[w->name] == -1) {
      pathdists[w->name] = pathdists[n->name] + 1;
      tmpqueue.push_back(w);
    }

    if (n->nchildren == 2) {
      Node *u = n->children[0];
      Node *v = n->children[1];

      if (pathdists[u->name] == -1) {
	pathdists[u->name] = pathdists[n->name] + 1;
	tmpqueue.push_back(u);
      }
      
      if (pathdists[v->name] == -1) {
	pathdists[v->name] = pathdists[n->name] + 1;
	tmpqueue.push_back(v);
      }
    }
  }

  sizequeuerevert=queue.size();
 
}

float SprNbrProposer::calcPropRatio(Tree *tree)
{

  float logratio=log(sizequeue)-log(sizequeuerevert);
  return logratio;
}



float NniProposer::calcPropRatio(Tree *tree)
{
  return 0;
}


void SprNbrProposer::pickNewSubtree()
{
    const Tree *tree = basetree;
    assert(basetree->nnodes >= 5);
    
    // find subtree (a) to cut off (any node that is not root or child of root)
    int choice;
    do {
        choice = irand(tree->nnodes);
    } while (tree->nodes[choice]->parent == NULL ||
             tree->nodes[choice]->parent->parent == NULL);
    Node *a = tree->nodes[choice];
    subtree = a;
    
    // find sibling (b) of a
    Node *c = a->parent;
    const int bi = (c->children[0] == a) ? 1 : 0;
    Node *b = c->children[bi];
    
    // uninitialize path distances
    pathdists.clear();
    for (int i=0; i<tree->nnodes; i++)
      pathdists.push_back(-1);
    
    // setup path distances and queue
    pathdists[a->name] = 0;
    pathdists[c->name] = 0;
    pathdists[b->name] = 0;
    queue.clear();
    list<Node*> tmpqueue;
    tmpqueue.push_back(c);
    tmpqueue.push_back(b);

    // traverse tree via depth first traversal
    while (tmpqueue.size() > 0) {
        Node *n = tmpqueue.front();
        tmpqueue.pop_front();
        
        // do not traverse beyond radius
        if (pathdists[n->name] >= radius)
            continue;

        // queue only valid new branch points:
        // n must not be root, a, descendant of a, c (parent of a), or  
        // b (sibling of a)
        if (n->parent && n != subtree && n != b && n != c) {
            queue.push_back(n);
        }

        // queue up unvisited neighboring edges
        Node *w = n->parent;

        if (w && pathdists[w->name] == -1) {
            pathdists[w->name] = pathdists[n->name] + 1;
            tmpqueue.push_back(w);
        }

        if (n->nchildren == 2) {
            Node *u = n->children[0];
            Node *v = n->children[1];

            if (pathdists[u->name] == -1) {
                pathdists[u->name] = pathdists[n->name] + 1;
                tmpqueue.push_back(u);
            }

            if (pathdists[v->name] == -1) {
                pathdists[v->name] = pathdists[n->name] + 1;
                tmpqueue.push_back(v);
            }
        }
    }
}

  ///////////////////////////LocalChangeProposer

SubtreeSlideProposer::SubtreeSlideProposer(int niter) :
    NniProposer(niter)
{
}


void SubtreeSlideProposer::propose(Tree *tree)
{
  performSubtreeSlide(tree, &m, &mstar);  
}


float SubtreeSlideProposer::calcPropRatio(Tree *tree)
{
  //printf("\non passe ds le Prop RAtio du subtreeSlide proposer\n");
  float logratio=3*log(mstar)-3*log(m);  
  return logratio;
  
}



  ///////////////////////////BranchLengthProposer, it changes only a little bite a random edge of the gene tree

BranchLengthProposer::BranchLengthProposer(int niter) :
    NniProposer(niter)
{
}

void BranchLengthProposer::propose(Tree *tree)
{

  performBranchLength(tree, &m, &mstar);  
}


float BranchLengthProposer::calcPropRatio(Tree *tree)
{
  float logratio=log(mstar)-log(m);  
  return logratio;


}



//=============================================================================
// Recon root proposer

void ReconRootProposer::propose(Tree *tree)
{
    const float rerootProb = 1.0;
    // propose new tree
    proposer->propose(tree);
    
    // reroot tree if stree is given
    if (frand() < rerootProb) {
        oldroot1 = tree->root->children[0];
        oldroot2 = tree->root->children[1];
        
        if (stree != NULL) {
            reconRoot(tree, stree, gene2species);
        }
    } else {
        oldroot1 = NULL;
        oldroot2 = NULL;
    }
}

void ReconRootProposer::revert(Tree *tree)
{
    // undo topology change
    if (oldroot1)
        tree->reroot(oldroot1, oldroot2);
    
    proposer->revert(tree);
}


//=============================================================================
// Dup/Loss proposer

DupLossProposer::DupLossProposer(TopologyProposer *proposer, 
                                 SpeciesTree *stree, int *gene2species,
                                 float dupprob, float lossprob,
                                 int quickiter, int niter) :
    proposer(proposer),
    quickiter(quickiter),
    niter(niter),
    iter(0),
    correctTree(NULL),
    correctSeen(false),
    stree(stree),
    gene2species(gene2species),
    dupprob(dupprob),
    lossprob(lossprob),
    recon(0),
    events(0),
    oldtop(NULL)
{
    doomtable = new double [stree->nnodes];
    calcDoomTable(stree, dupprob, lossprob, doomtable);

}


DupLossProposer::~DupLossProposer()
{
    delete [] doomtable;
}


void DupLossProposer::propose(Tree *tree)
{
    iter++;
    
    // do simple proposal if dup/loss probs are disabled
    if (dupprob < 0.0 || lossprob < 0.0 || quickiter <= 1) {
        proposer->propose(tree);
        return;
    }
    
    // save old topology
    oldtop = tree->copy();
       
    // recon tree to species tree
    recon.ensureSize(tree->nnodes);
    events.ensureSize(tree->nnodes);
    recon.setSize(tree->nnodes);
    events.setSize(tree->nnodes);
    
    ExtendArray<Tree*> trees(0, quickiter);
    ExtendArray<float> logls(0, quickiter);
    
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);
    double bestlogp = birthDeathTreePriorFull(tree, stree, recon, events, 
                                              dupprob, lossprob,
                                              doomtable,0.5);
    //use default value for the number of leaves at the root fix it later

    double sum = -INFINITY;

    // make many subproposals
    proposer->reset();
    for (int i=0; i<quickiter; i++) {
        proposer->propose(tree);
        // only allow unique proposals
        // but if I have been rejecting too much allow some non-uniques through
        if (uniques.has(tree) && trees.size() >= .1 * i) {
            proposer->revert(tree);
            continue;
        }

        reconcile(tree, stree, gene2species, recon);
        labelEvents(tree, recon, events);
        double logp = birthDeathTreePriorFull(tree, stree, recon, events, 
                                              dupprob, lossprob,
                                              doomtable,0.5);
	//use default value for the number of leaves at the root fix it later


        printLog(LOG_HIGH, "search: qiter %d %f %f\n", i, logp, bestlogp);
        
        Tree *tree2 = tree->copy();
        
        // save tree and logl
        trees.append(tree2);
        logls.append(logp);
        sum = logadd(sum, logp);
        
        if (logp > bestlogp)
            // make more proposals off this one
            bestlogp = logp;
        else
            proposer->revert(tree);
    }    
    
    // propose one of the subproposals 
    double choice = frand();
    double partsum = -INFINITY;
    
    for (int i=0; i<trees.size(); i++) {
        partsum = logadd(partsum, logls[i]);
        
        if (choice < exp(partsum - sum)) {
            // propose tree i
            printLog(LOG_MEDIUM, "search: choose %d %f %f\n", i, 
                     logls[i], exp(logls[i] - sum));
            tree->setTopology(trees[i]);            
            break;
        }
        
    }

    // add tree to unqiues
    uniques.insert(tree);

    // clean up subproposals
    for (int i=0; i<trees.size(); i++)
        delete trees[i];

}

void DupLossProposer::revert(Tree *tree)
{
    // do simple proposal if dup/loss probs are disabled
    if (dupprob < 0.0 || lossprob < 0.0 || quickiter <= 1) {
        proposer->revert(tree);
        return;
    }
    
    //printf("set oldtop\n");
    tree->setTopology(oldtop);
    delete oldtop;
}



//=============================================================================

UniqueProposer::~UniqueProposer()
{
    seenTrees.clear();
}


void UniqueProposer::propose(Tree *tree)
{
    iter++;
    for (int i=0;; i++) {
        printLog(LOG_HIGH, "search: unique trees seen %d (tries %d)\n", 
                 seenTrees.size(), i+1);

        proposer->propose(tree);
        
        if (seenTrees.insert(tree)) {
            // return new tree
            break;
        } else {
            if (i < ntries) {
                // revert and loop again	   
	      proposer->revert(tree);
            } else {
                // give up and return tree
                break;
            }
        }
    }
}


//=============================================================================
// get initial tree

// propose initial tree by Neighbor Joining
Tree *getInitialTree(string *genes, int nseqs, int seqlen, char **seqs)
{
    int nnodes = nseqs * 2 - 1;

    ExtendArray<int> ptree(nnodes);
    ExtendArray<float> dists(nnodes);
    Matrix<float> distmat(nseqs, nseqs);

    calcDistMatrix(nseqs, seqlen, seqs, distmat.getMatrix());
    neighborjoin(nseqs, distmat.getMatrix(), ptree, dists);

    Tree *tree = new Tree(nnodes);
    ptree2tree(nnodes, ptree, tree);
    tree->setLeafNames(genes);

    // special case
    if (nseqs == 2) {
        tree->nodes[0]->dist = dists[0];
        tree->nodes[1]->dist = dists[1];
    }

    return tree;
}


// propose initial tree and root by species
Tree *getInitialTree(string *genes, int nseqs, int seqlen, char **seqs,
                     SpeciesTree *stree, int *gene2species)
{
    Tree *tree = getInitialTree(genes, nseqs, seqlen, seqs);
    reconRoot(tree, stree, gene2species);

    return tree;
}


//=============================================================================
// search logging


void printSearchStatus(Tree *tree, SpeciesTree *stree, int *gene2species, int *recon, int *events, bool keepDupLoss, FILE *infoduploss, FILE *fileduplosslasttree,bool final)
{
  if (stree) {
    bool cleanupRecon = false;
    if (!recon) {
      recon = new int [tree->nnodes];
      cleanupRecon = true;
    }
    
    bool cleanupEvents = false;    
    if (!events) {
      events = new int [tree->nnodes];
      cleanupEvents = true;
    }    
	
    //   for (int i=0; i<tree->nnodes; i++) {
    // Node *node = tree->nodes[i];
      //printf("%d\t%s\n", node->name,(node->longname).c_str());
      //printf("its reconciliation is%d\n",recon[i]);
    //}

    int lossWGD=0;
    int losses = countLoss(tree, stree, recon, events, &lossWGD);
    printLog(LOG_LOW, "search: loss = %d\n", losses);
    printLog(LOG_LOW, "search: loss at the WGD = %d\n", lossWGD);

    // count dups
    int dups = 0;
    int dupWGD=0;
    for (int i=0; i<tree->nnodes; i++){
      if (events[i] == EVENT_DUP){
	dups++;
	if  (stree->nodes[recon[i]]->longname=="WGD_at"){
	  dupWGD++;
	}}
    }

    printLog(LOG_LOW, "search: dups = %d\n", dups);
    printLog(LOG_LOW, "search: dups at the WGD = %d\n", dupWGD);
    
    if (keepDupLoss){
    //we keep in a file the losses dupolications ...
    fprintf(infoduploss,"%d",losses);
    fprintf(infoduploss,"%d",lossWGD);
    fprintf(infoduploss,"%d",dups);
    fprintf(infoduploss,"%d",dupWGD);
    fprintf(infoduploss,"\n");
    }


    if (final){
      //we are at he last step ,ie last tree
      //we print in a file the losses duplication ......  only for the last tree
      fprintf(fileduplosslasttree,"%d",losses);
      fprintf(fileduplosslasttree,"%d",lossWGD);
      fprintf(fileduplosslasttree,"%d",dups);
      fprintf(fileduplosslasttree,"%d",dupWGD);

    }

        //assert(tree->assertTree());        
        if (cleanupRecon)
            delete [] recon;
        if (cleanupEvents)
            delete [] events;        
    }                
    
    if (isLogLevel(LOG_LOW))
    displayTree(tree, getLogFile());

}



void printLogTree(int loglevel, Tree *tree)
{
    if (isLogLevel(loglevel)) {
        printLog(loglevel, "tree: ");
        writeNewickTree(getLogFile(), tree, 0, true);
        printLog(loglevel, "\n");
    }
}


//////////////////////////////////////////////////////


void printTreeSampled(bool keepTreeSampled, FILE *filetrees, Tree *tree)
{
  if (keepTreeSampled){
     writeNewickTree(filetrees, tree, 0, true);
    fprintf(filetrees,"\n");
  }

}



///////////////////////////////////////////////////////






class Prob
{
public:
  double seqlk;
  double branchp;
  double topp;
  double logProbNotExtinct;
  double logp;


  double calcJoint(SpimapModel *model, Tree *tree)
  {
 	fflush(stdout);
        model->setTree(tree);
	seqlk = model->likelihood();
	branchp = model->branchPrior();
        topp = model->topologyPrior();
	logp = seqlk + branchp + topp - logProbNotExtinct ;
	fflush(stdout);

        return logp;
    }

  //  double calcJointWithoutTopp(SpimapModel *model, Tree *tree, float logProbNotExtinct)
  
  double calcJointWithoutTopp(SpimapModel *model, Tree *tree)
  { //we don t compute the topology prior
  
    model->setTree(tree);
    seqlk = model->likelihood();
    branchp = model->branchPrior();
    logp = seqlk + branchp + topp -  logProbNotExtinct;  

    return logp;
  }


 
 
  double calcJointWithBranchOptimization(SpimapModel *model, Tree *tree)
  {  
    //we don t compute the topology prior
    //we just optimize the branch length using original algorithm of Matt
    model->setTree(tree);
    seqlk = model->likelihoodWithOptimization();
    branchp = model->branchPrior();
    logp = seqlk + branchp + topp -  logProbNotExtinct;    

    return logp;
  }




};


void printLogProb(int loglevel, Prob *prob)
{
    printLog(loglevel, "search: lnl    = %f\n", prob->logp);
    printLog(loglevel, "search: seqlk  = %f\n", prob->seqlk);
    printLog(loglevel, "search: branch = %f\n", prob->branchp);          
    printLog(loglevel, "search: top    = %f\n", prob->topp);
    printLog(loglevel, "search: logProbNotExtinct = %f\n", prob->logProbNotExtinct);
}


TreeSearchClimb::TreeSearchClimb(SpimapModel *model, MixProposer *proposer, MixProposer *proposer2 ) :
    model(model),
    proposer(proposer),
    proposer2(proposer2)
{
}


TreeSearchClimb::~TreeSearchClimb()
{}


Tree *TreeSearchClimb::search(Tree *initTree, string *genes, 
			      int nseqs, int seqlen, char **seqs, string outputprefix, int method, bool keepTreeSampled, bool keepDupLoss)
{
    Tree *toptree = NULL;
    double logp = -INFINITY, nextlogp, logpuseless;
    Tree *tree = NULL;
    Timer correctTimer;    
    Timer proposalTimer;

    Prob prob;
    bool accept;
    double logPropRatio, logProbNotExtinct, logDoomedAtRoot;
    //    double q=0.5;//be careful it has to be the same q as in birthTreePrior2, fix it later

    double seqlk =0, branchp =0, topp =0, nextseqlk, nextbranchp, nexttopp;
    bool final=false;

    SpimapModel *model=getmodel();
    float q=model->getq();


    logDoomedAtRoot=model->getdoomtable()[model->getSpeciesTree()->root->name];
 
    logProbNotExtinct = log(1-exp(logDoomedAtRoot)) -  log(1-(1-q)*exp(logDoomedAtRoot)) ;
    prob.logProbNotExtinct= logProbNotExtinct;

    // setup search debug: testing against know correct tree
    Tree *correct = proposer->getCorrect();
    double correctLogp = -INFINITY;

    if (correct) {
        // determine probability of correct tree
        parsimony(correct, nseqs, seqs); // get initial branch lengths
        correctLogp = prob.calcJoint(model, correct);
        printLog(LOG_LOW, "search: correct tree lnl = %f\n", correctLogp);
    }

    
    // determine initial tree topology
    if (initTree == NULL)
	 tree = getInitialTree(genes, nseqs, seqlen, seqs,
                            model->getSpeciesTree(), 
                            model->getGene2species());


    if (initTree != NULL)
      tree=initTree->copy();


    // special cases (1 and 2 leaves)
    if (nseqs < 3) {
        return tree->copy();
    }
    

    // calc probability of initial tree
    parsimony(tree, nseqs, seqs); // get initial branch lengths
    logp = prob.calcJoint(model, tree);
    seqlk=prob.seqlk;
    branchp=prob.branchp;
    topp=prob.topp;

    toptree = tree->copy();



    // log initial tree
    printLog(LOG_LOW, "search: initial\n");
    printLogProb(LOG_LOW, &prob);
    printLogTree(LOG_LOW, tree);


    FILE *filetrees=NULL;
    if (keepTreeSampled){
      string outTreeSampledFile = outputprefix  + ".treesampled";
      filetrees=fopen(outTreeSampledFile.c_str(), "w");
    }else{
      filetrees=NULL;
    }

    printTreeSampled(keepTreeSampled, filetrees, tree);


    FILE *fileduploss=NULL;
    if (keepDupLoss){
      //we print all the history of the losses , duplication...corresponding to the different trees 
      string outDupLossFile = outputprefix  + ".duploss";
      fileduploss=fopen(outDupLossFile.c_str(), "w");
    }else{
      fileduploss=NULL;
    }


    //we also  print in a file the losses, duplications... of the last tree
    string duplosslasttreeFile = outputprefix  + ".duplosslasttree";
    FILE *fileduplosslasttreeFile=fopen(duplosslasttreeFile.c_str(), "w");  

    printSearchStatus(tree, model->getSpeciesTree(), model->getGene2species(), model->recon, model->events,keepDupLoss,fileduploss, fileduplosslasttreeFile,final);
    

    int naccept = 0;
    int nreject = 0;
    // search loop
    proposer->reset();
    
    fflush(stdout);
   
    for (int i=0; proposer->more(); i++) {


      printf("iteration  i of propose%d\n",i);

      printLog(LOG_LOW, "first stage :search iter %d\n", i);
    
      // propose new tree 
      proposalTimer.start();
      proposer->propose(tree);
	
      proposer->testCorrect(tree);
      proposal_runtime += proposalTimer.time();
        
      // calculate probability of proposal
      nextlogp = prob.calcJoint(model, tree);
      nextseqlk=prob.seqlk;
      nextbranchp=prob.branchp;
      nexttopp=prob.topp;

      accept=0;
      if (method==1){
	//MCMC
	logPropRatio=proposer->calcRatio(tree);
	accept = ((nextlogp > logp) ||  (frand()<exp(nextlogp-logp+logPropRatio)));
      }else{
	  //MAP ie maximum a posteriori
	accept = (nextlogp > logp);
      }
      

      // log proposal
      if (accept)
	printLog(LOG_LOW, "search: accept\n");
      else
	printLog(LOG_LOW, "search: reject\n");
      printLogProb(LOG_LOW, &prob);
      printLogTree(LOG_LOW, tree);

      // log for search debug
      if (correct && (nextlogp >= correctLogp || 
		      tree->sameTopology(correct)))
        {
	  printLog(LOG_LOW, "search: correct tree time = %f\n", 
		   correctTimer.time());
	  printLog(LOG_LOW, "search: correct tree logl = %f\n", 
		   correctLogp);
	  printLog(LOG_LOW, "search: correct tree is best = %d\n",
		   int(tree->sameTopology(correct)));
	  printLog(LOG_LOW, "search: correct tree better logl = %f\n", 
		   nextlogp);
	  
	  correct = NULL;
        }



      // act on acceptance
      if (accept) {
	naccept++;
	proposer->accept(true);
	logp = nextlogp;
	seqlk=nextseqlk;
	branchp=nextbranchp;
	topp=nexttopp;
	    
	delete toptree;
	toptree = tree->copy();
	      	    
	printSearchStatus(tree, model->getSpeciesTree(), model->getGene2species(), model->recon, model->events,keepDupLoss,fileduploss,fileduplosslasttreeFile,final);	    
	printTreeSampled(keepTreeSampled,filetrees,tree);	   

      } else {           
	// display rejected tree
	if (isLogLevel(LOG_MEDIUM))
	  printSearchStatus(tree, model->getSpeciesTree(), model->getGene2species(), model->recon, model->events,keepDupLoss,fileduploss,fileduplosslasttreeFile, final);
            	  
	nreject++;
	proposer->accept(false); 	     	     

	delete tree;
	tree=toptree->copy();

	// reject, undo topology change 
	printTreeSampled(keepTreeSampled,filetrees,tree);
        
      }


      printLog(LOG_LOW, "\n");


      /// SECOND STAGE or THIRD STAGE now at random

      prob.topp=topp;

      //printf("\nwe begin the Branchlength change stage\n");

      if (frand()<0.2){
	//SECOND STAGE
	//we propose little changes on the branchs lengths

	for (int k=0; k<(tree->nnodes-1); k++) {
	  printLog(LOG_LOW, "second stage :search iter %d\n", k);
	  proposer2->propose(tree);	 
	  proposer2->testCorrect(tree);

	  //we already have the topology prior
	  nextlogp = prob.calcJointWithoutTopp(model, tree);
	  nextseqlk=prob.seqlk;
	  nextbranchp=prob.branchp;
	  
	  accept=0;
	  if (method==1){
	  //MCMC
	    logPropRatio=proposer2->calcRatio(tree);
	    accept = ((nextlogp > logp) ||  (frand()<exp(nextlogp-logp+logPropRatio)));
	  }else{
	  //MAP ie maximum a posteriori
	    accept = (nextlogp > logp);
	  }
	  
	  // log proposal
	  if (accept)
            printLog(LOG_LOW, "search: accept\n");
	  else
            printLog(LOG_LOW, "search: reject\n");
	  printLogProb(LOG_LOW, &prob);
	  printLogTree(LOG_LOW, tree);

	  // act on acceptance
	  if (accept) {
            naccept++;
	    logp = nextlogp;
	    seqlk=nextseqlk;
	    branchp=nextbranchp;
	    
            delete toptree;

	    toptree = tree->copy();	   
	    printSearchStatus(tree, model->getSpeciesTree(), model->getGene2species(), model->recon, model->events,keepDupLoss,fileduploss,fileduplosslasttreeFile,final);
	    
	  }else{	    

	    nreject++;	    	    
	    delete tree;
	    tree = toptree->copy();
	    
	  }

	}
	//end SECOND STAGE

      }else{

	/// THIRD STAGE 
	//we use the optimization of Matt for the branch length

	printLog(LOG_LOW, "third stage :search iter 0\n");
	proposer2->propose(tree);	 
	proposer2->testCorrect(tree);

	//we already have the topology prior
	nextlogp = prob.calcJointWithBranchOptimization(model, tree);
	nextseqlk=prob.seqlk;
	nextbranchp=prob.branchp;

	accept=0;
	if (method==1){
	  //MCMC
	  logPropRatio=proposer2->calcRatio(tree);	
	  accept = ((nextlogp > logp) ||  (frand()<exp(nextlogp-logp+logPropRatio)));
	}else{
	  //MAP ie maximum a posteriori
	  accept = (nextlogp > logp);
	}
	
	// log proposal
	if (accept)
	  printLog(LOG_LOW, "search: accept\n");
	else
	  printLog(LOG_LOW, "search: reject\n");
	printLogProb(LOG_LOW, &prob);
	printLogTree(LOG_LOW, tree);

	// act on acceptance
	if (accept) {
	  naccept++;
	  logp = nextlogp;
	  seqlk=nextseqlk;
	  branchp=nextbranchp;
	    
	  delete toptree;	  
	  toptree = tree->copy();
	   
	  printSearchStatus(tree, model->getSpeciesTree(), model->getGene2species(), model->recon, model->events,keepDupLoss,fileduploss,fileduplosslasttreeFile,final);
	    
	}else{	    

	  nreject++;	    
	  delete tree;
	  tree = toptree->copy();
	    
	}

	//end THIRD STAGE

      }

      printTreeSampled(keepTreeSampled,filetrees,tree);

    }

    ///////////////////////////////////////////////////
    
    // print final log messages
    printLog(LOG_LOW, "accept rate: %f\n", naccept / double(naccept+nreject));
    logpuseless=prob.calcJointWithoutTopp(model, tree);

    //be careful, we already had saved thebest logp and the corresponding seqlk branchp topp 
    //however we had to run again calcjoint above in order to update model->recon end model->events
    //which had not been saved 
    //otherwise we have model->recon and model->events not corresponding to the best tree 
    //so for the following we take the thing saved
    // and we take the new model->recon and model->events

    printLog(LOG_LOW, "search: final\n");

    prob.seqlk=seqlk;
    prob.branchp=branchp;
    prob.topp=topp;
    prob.logp=logp;
    printLogProb(LOG_LOW, &prob);
    printLog(LOG_LOW, "double check %f\n", logp);
    
    final=true;
    //in order to write in a file the losses and duplications corresponding to the last tree
  
    printSearchStatus(tree, model->getSpeciesTree(), model->getGene2species(), model->recon, model->events,keepDupLoss,fileduploss,fileduplosslasttreeFile,final);    

    fclose(fileduplosslasttreeFile);

    if (keepTreeSampled){
      fclose(filetrees);}

    if (keepDupLoss){
      fclose(fileduploss);}

    // print the full log likelihood in the file 
    string logpFile = outputprefix  + ".completeloglikelihood";
    FILE *filelogpFile=fopen(logpFile.c_str(), "w");      

    fprintf(filelogpFile,"%e",prob.logp);
    
    // print the probability of the topology  in the file
    string toppFile = outputprefix  + ".topologyprob";
    FILE *filetoppFile=fopen(toppFile.c_str(), "w");      

    fprintf(filetoppFile,"%e",prob.topp);

    fclose(filelogpFile);
    fclose(filetoppFile);
    
    delete tree;
    return toptree;

    }


/*

extern "C" {


// Calculate the likelihood of a tree
Tree *searchClimb(int niter, int quickiter,
		  int nseqs, char **gene_names, char **seqs,
		  int nsnodes, int *pstree, float *sdists,
		  int *gene2species,
		  float *sp_alpha, float *sp_beta, float generate,
		  float pretime_lambda, float birth, float death,
		  float gene_alpha, float gene_beta,
		  float *bgfreq, float kappa,
		  int nsamples, bool approx)
{
    int nnodes = 2*nseqs - 1;

    setLogLevel(LOG_MEDIUM);

    // create stree
    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();
    stree.setDists(sdists);
    
    // params
    SpidirParams params(nsnodes, NULL, sp_alpha, sp_beta, 
			gene_alpha, gene_beta, pretime_lambda);
    
    // make gene names
    string genes[nseqs];
    //char s[101];
    for (int i=0; i<nseqs; i++) {
	//snprintf(s, 100, "%d", i);
	genes[i] = gene_names[i];
    }

    // model
    const int maxiter = 2;
    int seqlen = strlen(seqs[0]);
    //  SpimapModel model(nnodes, &stree, &params, gene2species,
    //		      pretime_lambda, birth, death, nsamples, approx, false);
    
    SpimapModel model(nnodes, &stree, &stree, &params, gene2species,
		      pretime_lambda, birth, death, nsamples, approx, false);
    
    model.setLikelihoodFunc(new HkySeqLikelihood(nseqs, seqlen, seqs, 
                                                 bgfreq, kappa, maxiter));
    // proposers
    NniProposer nni(niter);
    SprProposer spr(niter);
    MixProposer mix(niter);
    mix.addProposer(&nni, .5);
    mix.addProposer(&spr, .5);
    ReconRootProposer rooted(&mix, &stree, gene2species);
    UniqueProposer unique(&rooted, niter);
    DupLossProposer proposer(&mix, &stree, gene2species, 
			     birth, death,
                             quickiter, niter);
    
    TreeSearchClimb search(&model, &proposer);

    Tree *tree = search.search(NULL, genes, nseqs, seqlen, seqs);
    
    return tree;
}

} // extern "C"

*/
} // namespace spidir
