/*=============================================================================

  Rabier Charles-Elie
  Copyright 2007-2011

  Newick tree reading/writing

=============================================================================*/

// c++ headers
#include <assert.h>
#include <stdio.h>

// spidir headers
#include "Tree.h"
#include "phylogeny.h"
#include "parsing.h"
#include "logging.h"
#include "WGD.h"
#include "model_params.h"


namespace spidir {


// function to give to all the nodes corresponding to the different
// WGD the largest indices
int nodeNameCmpWGD(const void *_a, const void *_b)
{
    Node *a = *((Node**) _a);
    Node *b = *((Node**) _b);

    if ((a->longname == "WGD_before") || (a->longname == "WGD_at")){
      if ((b->longname == "WGD_before") || (b->longname == "WGD_at"))
	return 0;
      else	
	return 1;}
    else{
      if ((b->longname == "WGD_before") || (b->longname == "WGD_at"))
	return -1;
      else	
	return 0;
    }
}



// function to remove WGDnodes in order to obtain a normal tree
//this function do two things
//first it gives the WGD nodes the largest node names
//then it return the new WGD tree
// and also the tree without WGD (modify argument tree)  
//this way, the tree with WGD and without WGD have exactly the same node names for the 
//non WGD nodes

SpeciesTree *removeWGDnodes(SpeciesTree *tree)
{
  //we give the biggest node names for the WGD nodes
  qsort((void*) tree->nodes.get(), tree->nodes.size(), 
	sizeof(Node*), nodeNameCmpWGD);
 

  // update names
  for (int i=0; i<tree->nnodes; i++)
    tree->nodes[i]->name = i;

  SpeciesTree * WGDtree = tree->copy();
 
  //we will consider the different WGD i
  for (int i=0; i<tree->nWGD; i++){
   
    //we need to obtain parent of WGD_before
    Node *parent = tree->nodes[tree->theWGD[i]->WGD_before->parent->name];

    //we need to obtain the node WGD_after of the WGD
    Node *after=tree->nodes[tree->theWGD[i]->WGD_after->name]; 

    //need to be careful if there are more than one WGD on the segment considered
    //index will be the index of the most recent WGD on the segment
    int index=0;

    double distWGD=tree->theWGD[0]->WGD_after->dist;
    //dist WGD is total length of the segment where al the the WGD takes place 

    if (i>0){
      index=i;
      distWGD=tree->theWGD[index]->WGD_after->dist;

      while  (tree->nodes[tree->theWGD[index-1]->WGD_before->parent->name]==
	      tree->nodes[tree->theWGD[index]->WGD_at->name]){
	--index;
	distWGD=distWGD + tree->theWGD[index]->WGD_after->dist;
	
	if (index==0)
	  break;
      }
      after=tree->nodes[tree->theWGD[index]->WGD_after->name];
      
    }
    
    distWGD=distWGD+tree->theWGD[i]->WGD_before->dist;
    
    //we search the index j corresponding to the children of parent
    int j = 0;
    for (; j<=parent->nchildren; j++) {
      if (parent->children[j] == tree->nodes[tree->theWGD[i]->WGD_before->name])
	break;
    }

    after->parent=parent;
    after->dist=distWGD;
    parent->children[j]=after;
  }

  // the WGD nodes have the biggest indices
  for (int i=tree->nnodes; i>tree->nnodes-2*tree->nWGD; i--) {  
    tree->nodes.pop();
  }

  tree->nnodes=tree->nodes.size();
  tree->nWGD=0;
  delete [] tree->theWGD;

  return WGDtree;
}

  ///////////////////////////////////////////////////////////////////////////////////

//function to prepare parameters for WGD
//assumes that in the WGDtree, all the nodes WGD_before and WGD_at
//have a name larger than the other nodes

//function which extend params which contained rates only for nodes which were not WGDnodes 
//if smik is a node with a WGD just above
//we keep the same rate (for WGD-at and WGD-before) as the one for smik read in the file .params

void extendRateParamToWGDnodes(SpidirParams *params, SpeciesTree *WGDstree)
{   
    params->extendSp_alphaSp_beta(WGDstree->nnodes);

    int i=0;
    for (; i<WGDstree->nWGD; i++) {

	params->sp_alpha[WGDstree->theWGD[i]->WGD_at->name]=params->sp_alpha[WGDstree->theWGD[i]->WGD_after->name];
	
     	params->sp_alpha[WGDstree->theWGD[i]->WGD_before->name]=params->sp_alpha[WGDstree->theWGD[i]->WGD_at->name];
	
	params->sp_beta[WGDstree->theWGD[i]->WGD_at->name]=params->sp_beta[WGDstree->theWGD[i]->WGD_after->name];
	
	params->sp_beta[WGDstree->theWGD[i]->WGD_before->name]=params->sp_beta[WGDstree->theWGD[i]->WGD_at->name];
	
      }
 
}


  //function which find the best location for a node, ie before at or after the WGD
  //criteria : maximize the number of duplications which happened at the WGD
int  WGDreconcile_onebranch(SpeciesTree *WGDstree,  int *recon_noWGD, int *events, Node *nodeWGDafter,Node *node,int *best)
{
  int dup_at; 

  if ((recon_noWGD[node->name]!=nodeWGDafter->name) || (events[node->name]!=EVENT_DUP)){

    // node is  a pseudoLeaf
    best[node->name]=AFTER;
    dup_at=0;
    //losses_at=1;
  }else
    {//node is not a pseudoleaf
      
      int dup_at0;
      int dup_at1;
      dup_at0=WGDreconcile_onebranch(WGDstree, recon_noWGD, events, nodeWGDafter ,node->children[0],best);
      
      dup_at1=WGDreconcile_onebranch(WGDstree, recon_noWGD, events, nodeWGDafter ,node->children[1],best);
      
      if ((dup_at0+dup_at1)>1) {
	  best[node->name]=BEFORE ;
	  dup_at=dup_at0+dup_at1;
	}else
	{
	  best[node->name]=AT ;
	  dup_at=1;
	}
      
    }

return dup_at;

}


  //function which uses the information best (obtained with  WGDreconcile_onebranch), and changes recon based on  best
  void WGDreconReset(SpeciesTree *WGDstree,  int *recon, int *recon_noWGD, int *events, Node *nodeWGDafter,Node *node,int *best)
{
  bool pseudoleaf;
  pseudoleaf=0;

  if ((recon_noWGD[node->name]!=nodeWGDafter->name) || (events[node->name]!=EVENT_DUP)){
      pseudoleaf=1;
    }

  if (!pseudoleaf){

    if (best[node->name]==AT){
      //reconciliation with the node WGDat
      recon[node->name]=nodeWGDafter->parent->name;
      // events[node->name]= EVENT_DUP;
    }

    if (best[node->name]==BEFORE){ 

      //reconciliation with the node WGDbefore
      recon[node->name]=nodeWGDafter->parent->parent->name;
      //events[node->name]= EVENT_DUP;
      WGDreconReset(WGDstree, recon, recon_noWGD, events, nodeWGDafter, node->children[0], best);
      WGDreconReset(WGDstree, recon, recon_noWGD, events, nodeWGDafter, node->children[1], best);
      
    }

  }


}


}// end of namespace spidir
