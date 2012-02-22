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
  writeNewickTree(stdout, tree,true);
    // undo topology change
  performNni(tree, nodea, nodeb); 
  writeNewickTree(stdout, tree,true);

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
// Mixuture of Proposers

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
   //print queue
    printf("\nthis is our non random queue\n");
    list<Node*>::iterator indice=queue.begin();
    for (int i=0; i<queue.size(); i++){
      printf("\nname of the node %d\n",(*indice)->name);
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

    printf("\nthis is our random queue\n");
    list<Node*>::iterator indices=randomqueue.begin();
    for (int i=0; i<randomqueue.size(); i++){
      printf("\nname of the node %d\n",(*indices)->name);
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

LocalChangeProposer::LocalChangeProposer(int niter) :
    NniProposer(niter)
{
}


void LocalChangeProposer::propose(Tree *tree)
{
  performLocalChange(tree, &m, &mstar);  
}


float LocalChangeProposer::calcPropRatio(Tree *tree)
{
  float result=(mstar/m)*(mstar/m)*(mstar/m);  
  return result;
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
                                              doomtable);


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
                                              doomtable);
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
	      printf("on effectue un revert ds Uniqueproposer");
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


void printSearchStatus(Tree *tree, SpeciesTree *stree, int *gene2species, int *recon, int *events,FILE *infoduploss)
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
	
    for (int i=0; i<tree->nnodes; i++) {
      Node *node = tree->nodes[i];
      printf("%d\t%s\n", node->name,(node->longname).c_str());
      printf("its reconciliation is%d\n",recon[i]);
    }

    int lossWGD=0;
    int losses = countLoss(tree, stree, recon, events, &lossWGD);
    printf("we have found %d losses \n",losses);
    printLog(LOG_LOW, "search: loss = %d\n", losses);
    printf("we have found %d losses at the WGD \n",lossWGD);
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

    printf("we have found %d dup \n",dups);
    printf("we have found %d dup at the WGD \n",dupWGD);
    printLog(LOG_LOW, "search: dups = %d\n", dups);
    printLog(LOG_LOW, "search: dups at the WGD = %d\n", dupWGD);
    
    fprintf(infoduploss,"%d",losses);
    fprintf(infoduploss,"%d",lossWGD);
    fprintf(infoduploss,"%d",dups);
    fprintf(infoduploss,"%d",dupWGD);
    fprintf(infoduploss,"\n");

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


class Prob
{
public:
    double seqlk;
    double branchp;
    double topp;
    double logp;

  double calcJoint(Model *model, Tree *tree, double *theseqlik, double *thelogp)
  {
        printf("now starting setTree\n");
	fflush(stdout);

        model->setTree(tree);

	seqlk = model->likelihood();
	printf("\noh voici seqlklikelihood %f\n",seqlk);
	*theseqlik=seqlk;
	printf("\nre voici expseqlklikelihood %f\n",*theseqlik);

	branchp = model->branchPrior();

	printf("Branch Prior finished");
	fflush(stdout);

	///affichea virer
	for (int i=0; i<tree->nnodes; i++) {
	  Node *node = tree->nodes[i];
	  printf("%d\t%s\n", node->name,(node->longname).c_str());
	  for (int i=0; i<node->nchildren; i++){
	    printf("\t%d\t%s\n", node->children[i]->name, (node->children[i]->longname).c_str());
	  }
	}
   //endafficheavirer

        topp = model->topologyPrior();
        logp = seqlk + branchp + topp;
	*thelogp=logp;
	printf("sum finished");
	fflush(stdout);

        return logp;
    }
};


void printLogProb(int loglevel, Prob *prob)
{
    printLog(loglevel, "search: lnl    = %f\n", prob->logp);
    printLog(loglevel, "search: seqlk  = %f\n", prob->seqlk);
    printLog(loglevel, "search: branch = %f\n", prob->branchp);          
    printLog(loglevel, "search: top    = %f\n", prob->topp);
}



//============================================================


/*


void calcSumForHarmonic(double likseq, double *harmonic1)
{   

*harmonic1= *harmonic1 + 1/likseq;

}


void calcFinalHarmonic(int niter, double *harmonic1)
{
  
  //before
  // *harmonic1=1/(*harmonic1/niter);
  //end before

  double result;

  result= -log(90) + log(1/(*harmonic1/niter)) ;
  
  *harmonic1=exp(result);


}


void calcSumForHarmonic4(double likseq, double *sumnumerator, double *sumdenominator)
{   

  double value=0;
  double delta=0;


  for (int i=1; i<101; ++i){    
    value += 0.01;
    sumnumerator[i] +=  likseq/(value*delta + (1-delta)*likseq); 
    sumdenominator[i] += 1/(value*delta + (1-delta)*likseq);

    //    printf("\nvoici i %d\n",i);
    // printf("\nvoici sumnumerator %f\n",sumnumerator[10]);
    //printf("\n o lala voici sumdenominator %f\n",sumdenominator[i]);
  


  }


}


void calcFinalHarmonic4(int niter, double *harmonic4, double *sumnumerator, 
double *sumdenominator)
{

  double value=0;
  double delta=0.1;
  double num, denom, diff, olddiff;
  *harmonic4=3;
  olddiff=5;


  for (int i=1; i<101; ++i){ 

    //  printf("\nvoici i %d\n",i);
    //printf("\nvoici sumnumerator %f\n",sumnumerator[i]);
    //printf("\nvoici sumdenominator %f\n",sumdenominator[i]);

    value += 0.01;
    num=delta*niter/(1-delta) + sumnumerator[i]; 
    //    denom= delta*niter/(value*(1-delta)) + sumdenominator[i];
    denom= value*delta*niter/(1-delta) + sumdenominator[i];
    diff= (value-num/denom)*(value-num/denom);

    // printf("\nvoici diff %f\n",diff);

    if ( min(diff,olddiff)!= olddiff ){
      *harmonic4=value;
      olddiff=diff;
    }

  }


}


*/







//=============================================================================
// Search Climb

//TreeSearchClimb::TreeSearchClimb(Model *model, TopologyProposer *proposer) :
  TreeSearchClimb::TreeSearchClimb(Model *model, MixProposer *proposer) :
    model(model),
    proposer(proposer)
{
}


TreeSearchClimb::~TreeSearchClimb()
{}


Tree *TreeSearchClimb::search(Tree *initTree, string *genes, 
			      int nseqs, int seqlen, char **seqs, string outputprefix,int method)
{
    Tree *toptree = NULL;
    double toplogp = -INFINITY, nextlogp;
    Tree *tree = initTree;
    Timer correctTimer;    
    Timer proposalTimer;

    Prob prob;
    double theseqlik, theoldseqlik, thelogp, theoldlogp;

    
    // setup search debug: testing against know correct tree
    Tree *correct = proposer->getCorrect();
    double correctLogp = -INFINITY;
    if (correct) {
        // determine probability of correct tree
        parsimony(correct, nseqs, seqs); // get initial branch lengths
        correctLogp = prob.calcJoint(model, correct, &theseqlik, &thelogp);
        printLog(LOG_LOW, "search: correct tree lnl = %f\n", correctLogp);
    }


    
    // determine initial tree topology
    //   if (initTree == NULL)
    //  tree = getInitialTree(genes, nseqs, seqlen, seqs,
    //                        model->getSpeciesTree(), 
    //                        model->getGene2species());


    // special cases (1 and 2 leaves)
    if (nseqs < 3) {
        return tree->copy();
    }
    

    // calc probability of initial tree
    parsimony(tree, nseqs, seqs); // get initial branch lengths
    toplogp = prob.calcJoint(model, tree, &theseqlik, &thelogp);
    theoldseqlik=theseqlik;
    theoldlogp=thelogp;

    toptree = tree->copy();
    
    // log initial tree
    printLog(LOG_LOW, "search: initial\n");
    printLogProb(LOG_LOW, &prob);
    printLogTree(LOG_LOW, tree);

    string outTreeSampledFile = outputprefix  + ".treesampled";
    FILE *filetrees=fopen(outTreeSampledFile.c_str(), "w");
    writeNewickTree(filetrees, tree, 0, true);
    fprintf(filetrees,"\n");
 

    string outDupLossFile = outputprefix  + ".duploss";
    FILE *fileduploss=fopen(outDupLossFile.c_str(), "w");

    printSearchStatus(tree, model->getSpeciesTree(), model->getGene2species(), model->recon, model->events,fileduploss);
   
    int naccept = 0;
    int nreject = 0;
    // search loop
    proposer->reset();
    
    printf("\non est avt la boucle for\n");
    fflush(stdout);
   
    for (int i=0; proposer->more(); i++) {
        printLog(LOG_LOW, "search: iter %d\n", i);
    
        // propose new tree 
        proposalTimer.start();
	printf("\non est avt le propose de treesearchclimb\n");

	//avirer apres
	printf("\noh on est juste avt de changer\n");
	//	writeNewickTree(stdout, tree,true);
   

   printf("\n on affiche notre gene tree\n");
    for (int k=0; k<tree->nnodes; k++) {
      Node *node = tree->nodes[k];
      printf("%d\t%s\n", node->name,(node->longname).c_str());
      for (int j=0; j<node->nchildren; j++){
	printf("\t%d\t%s\n", node->children[j]->name, (node->children[j]->longname).c_str());
      }
    }
   writeNewickTree(stdout, tree,true);
	//end eavirer

        proposer->propose(tree);
	printf("\non est apres le propose de treesearchclimb\n");

	//avirer apres
  printf("\n on affiche notre gene tree apre le propose\n");
    for (int k=0; k<tree->nnodes; k++) {
      Node *node = tree->nodes[k];
      printf("%d\t%s\n", node->name,(node->longname).c_str());
      for (int j=0; j<node->nchildren; j++){
	printf("\t%d\t%s\n", node->children[j]->name, (node->children[j]->longname).c_str());
      }
    }
 writeNewickTree(stdout, tree,true);
 //end a virer apres
	
	printf("\nle voici notre i%d\n",i);
	fflush(stdout);

	//avirer apres
	printf("\noh on est juste apres avoir changer\n");
	//	writeNewickTree(stdout, tree,true);
	//end eavirer


	proposer->testCorrect(tree);
        proposal_runtime += proposalTimer.time();
        
        // calculate probability of proposal
        nextlogp = prob.calcJoint(model, tree, &theseqlik, &thelogp);

	//avirer apres
	printf("\noh on est juste apres avoir calcule le joint\n");
	fflush(stdout);
	//	writeNewickTree(stdout, tree,true);
	//end eavirer

	bool accept=0;
	if (method==1){
	  //MCMC
	float logPropRatio=proposer->calcRatio(tree);
	printf("\nvoici le logPropRatio%f\n", logPropRatio);
	accept = ((nextlogp > toplogp) ||  (frand()<exp(nextlogp-toplogp+logPropRatio)));
	}else{
	  //ML ie maximum a posteriori
	 accept = (nextlogp > toplogp);
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
            toplogp = nextlogp;
            delete toptree;
            toptree = tree->copy();
	   
	    //try something
	    if (theseqlik>theoldseqlik){
	      theoldseqlik=theseqlik;}
	    //end try

	    //try something
	    if (thelogp>theoldlogp){
	      theoldlogp=thelogp;}
	    //end try


	    printSearchStatus(tree, model->getSpeciesTree(), model->getGene2species(), model->recon, model->events,fileduploss);
	    
	    writeNewickTree(filetrees, tree, 0, true);
	    fprintf(filetrees,"\n");

        } else {           
            // display rejected tree
            if (isLogLevel(LOG_MEDIUM))
	      printSearchStatus(tree, model->getSpeciesTree(), model->getGene2species(), model->recon, model->events,fileduploss);
            	  
	      nreject++;
	      proposer->accept(false); 
 
	      //try something
	      if (theseqlik>theoldseqlik){
	      theoldseqlik=theseqlik;}
	      //end try

	      //try something
	      if (thelogp>theoldlogp){
	      theoldlogp=thelogp;}
	      //end try

	      delete tree;
	      tree=toptree->copy();

	      // reject, undo topology change
	      writeNewickTree(filetrees, tree, 0, true);
	      fprintf(filetrees,"\n");
        
        }

        printLog(LOG_LOW, "\n");
    }
    
    // print final log messages
    printLog(LOG_LOW, "accept rate: %f\n", naccept / double(naccept+nreject));

    // call probability break down again of top tree
    prob.calcJoint(model, toptree, &theseqlik, &thelogp);
    
    // log final tree
    printLog(LOG_LOW, "search: final\n");
    printLogProb(LOG_LOW, &prob);

    printSearchStatus(toptree, model->getSpeciesTree(), model->getGene2species(), model->recon, model->events,fileduploss);


    fclose(filetrees);
    fclose(fileduploss);

    // print the last likelihood of the sequence

    string loglikelihoodseqFile = outputprefix  + ".loglikelihoodseq";
    FILE *fileloglikelihoodseq=fopen(loglikelihoodseqFile.c_str(), "w");     
    fprintf(fileloglikelihoodseq,"%e",theoldseqlik);
    fclose(fileloglikelihoodseq);
    //

    // print the last likelihood of the sequence
    string logpFile = outputprefix  + ".completeloglikelihood";
    FILE *filelogpFile=fopen(logpFile.c_str(), "w");     
    fprintf(filelogpFile,"%e",theoldlogp);
    fclose(filelogpFile);
    //



    // clean up
    if (initTree == NULL)
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
