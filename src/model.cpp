/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gene tree evolution model

=============================================================================*/


#include "common.h"
#include "branch_prior.h"
#include "logging.h"
#include "Matrix.h"
#include "model.h"
#include "phylogeny.h"
#include "seq_likelihood.h"
#include "top_prior.h"
#include "WGD.h"
#include "treevis.h"

namespace spidir {



SpimapModel::SpimapModel(
    int nnodes, SpeciesTree *stree, SpeciesTree *stree_small,
    SpidirParams *params, 
    int *gene2species,
    float predupprob, float dupprob, float lossprob, 
    int nsamples, bool approx, bool useBranchPrior, float q) :
    
    Model(nnodes),
    nnodes(nnodes),
    stree(stree),
    stree_noWGD(stree_small),
    params(params),
    gene2species(gene2species),
    predupprob(predupprob),
    dupprob(dupprob),
    lossprob(lossprob),
    nsamples(nsamples),
    approx(approx),
    useBranchPrior(useBranchPrior),
    q(q)
{
    doomtable = new double [stree->nnodes]; 
    calcDoomTable(stree, dupprob, lossprob, doomtable);    

}


SpimapModel::~SpimapModel()
{
    delete [] doomtable;

    if (likelihoodFunc)
        delete likelihoodFunc;
}

void SpimapModel::setTree(Tree *_tree)
{
    tree = _tree;
    spidir::reconcile(tree, stree_noWGD, gene2species, recon_noWGD);
    labelEvents(tree, recon_noWGD, events);

    for (int j=0; j<tree->nnodes; j++) {
       recon[j]=recon_noWGD[j];     
     }

     WGDreconcile(tree->root,-1);    

}



  void SpimapModel::WGDreconcile(Node *node,int thelastWGD) 
{

  int lastWGD;
  lastWGD=thelastWGD;
  int dupat;

  if (recon_noWGD[node->name]!=thelastWGD){

    lastWGD=-1;

    for (int j=0; j<stree->nWGD; j++) {

      if (recon_noWGD[node->name]==stree->theWGD[j]->WGD_after->name){
	
	dupat=WGDreconcile_onebranch(stree, recon_noWGD, events, stree->theWGD[j]->WGD_after, node, best);
      	WGDreconReset(stree, recon, recon_noWGD, events, stree->theWGD[j]->WGD_after, node, best);
	lastWGD= stree->theWGD[j]->WGD_after->name;

      }

    }

  }

  for (int k=0; k<node->nchildren; k++) {
    WGDreconcile(node->children[k], lastWGD);
  }

  
}


  /*
  /////////////////////////////////////try to propose a new reconciliation for the WGD (ie note the most parsimonious)
//call the function WGDReconcileBottom(tree->root) 



  void SpimapModel::WGDReconcileBottom(Node *node) 
{


   for (int j=0; j<stree->nWGD; j++) {

     if ((recon[node->name] == stree->theWGD[j]->WGD_before->name)|| (recon[node->name] == stree->theWGD[j]->WGD_at->name))
       {
	recon[node->name]= stree->theWGD[j]->WGD_after->name;
       }

   }
   

  for (int k=0; k<node->nchildren; k++) {
  WGDReconcileBottom(node->children[k]);
  }


}


  void SpimapModel::WGDReconcileMiddle(Node *node) 
{

 for (int j=0; j<stree->nWGD; j++) {

   if ((recon[node->name] == stree->theWGD[j]->WGD_after->name) && (events[node->name] ==EVENT_DUP)  )
       {
	 recon[node->name]= stree->theWGD[j]->WGD_at->name;
       }
 }


 for (int k=0; k<node->nchildren; k++) {
   WGDReconcileBottom(node->children[k]);
 }


}

  getSpecSubtree(node, schild, recon, events, subleaves);



  ///////end try new reconciliation

  */



double SpimapModel::likelihood()
{
    double logp = 0.0;
    if (likelihoodFunc) {
        Timer timer;
        logp = likelihoodFunc->findLengths(tree);
        seq_runtime += timer.time();
    }
    return logp;
}



double SpimapModel::likelihoodWithOptimization()
{
    double logp = 0.0;
    if (likelihoodFunc) {
        Timer timer;
        logp = likelihoodFunc->findLengthsWithOptimization(tree);
        seq_runtime += timer.time();
    }
    return logp;
}




double SpimapModel::branchPrior()
{
    if (isNullParams(params) || !useBranchPrior) {
        return 0.0;
    }

    Timer timer;
    const float generate = -99; // integrate over gene rate
    double logp = spidir::branchPrior(tree, stree,
                                      recon, events, params,
                                      generate, predupprob, dupprob, lossprob,
                                      nsamples, approx);
    branch_runtime += timer.time();
    return logp;
}


double SpimapModel::topologyPrior()
{
    Timer timer;

    //printf("just before Treepriorfull");
	fflush(stdout);
    double logp = birthDeathTreePriorFull(tree, stree, recon, events, 
					  dupprob, lossprob, doomtable,q);

    top_runtime += timer.time();
    return logp;
}


//=============================================================================
// HKY sequence likelihood

HkySeqLikelihood::HkySeqLikelihood(int nseqs, int seqlen, char **seqs, 
                                   float *bgfreq, float tsvratio, int maxiter,
                                   double minlen, double maxlen) :
    nseqs(nseqs),
    seqlen(seqlen),
    seqs(seqs),
    bgfreq(bgfreq),
    tsvratio(tsvratio),
    maxiter(maxiter),
    minlen(minlen),
    maxlen(maxlen)
{}


double HkySeqLikelihood::findLengths(Tree *tree)
{ 


  //we only want the likelihood of the sequences given our tree 
  //no ML to find the branch lengths
     return  calcSeqProbHky(tree,nseqs,seqs, 
     		       bgfreq,tsvratio);



}



double HkySeqLikelihood::findLengthsWithOptimization(Tree *tree)
{ 


     return findMLBranchLengthsHky(tree, nseqs, seqs, bgfreq, 
  				  tsvratio, maxiter, minlen, maxlen);



}




} // namespace spidir

