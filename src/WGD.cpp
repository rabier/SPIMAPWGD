/*=============================================================================

  Matt Rasmussen # fixit
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
SpeciesTree *removeWGDnodes(SpeciesTree *tree)
{
  //we give the biggest node names for the WGD nodes
  qsort((void*) tree->nodes.get(), tree->nodes.size(), 
	sizeof(Node*), nodeNameCmpWGD);
 

 // update names
  for (int i=0; i<tree->nnodes; i++)
    tree->nodes[i]->name = i;


    
  //Tree * WGDtree = tree->copy(); // does not copy the WGD paramters
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






//function to prepare parameters for WGD
  SpidirParams  *paramforWGD(SpidirParams *params, SpeciesTree *WGDstree)
{


  //const int MAX_NAME = 51;
  // float param1, param2;
 float alpha = -1, beta = -1;
 ExtendArray<float> mu(0, 40);
 ExtendArray<float> sigma(0, 40);
 ExtendArray<string> names(0, 40);


 //parameters for the gene
 alpha = params->gene_alpha;
 beta = params->gene_beta;
 
 //parameters for  WGDstree

 for (int i=0; i<params->nsnodes; i++) {
 
   names.append(params->names[i]);
   mu.append(params->sp_alpha[i]);
   sigma.append(params->sp_beta[i]);

 }


 for (int i=params->nsnodes; i<WGDstree->nnodes; i++){
   names.append("");
   mu.append(0);
   sigma.append(0);

 }

  /*
 printf("voila params->nsnodes %d\n", params->nsnodes);
 printf("voila names.len %d\n", names.size());

 printf("voila WGDstree->nWGD %d\n",WGDstree->nWGD);
 //names.setCapacity(params->nsnodes+2*WGDstree->nWGD);
names.setCapacity(33);
mu.setCapacity(33);
sigma.setCapacity(33);


printf("voila names.len %d\n", names.size());
//mu.setCapacity(params->nsnodes+2*WGDstree->nWGD);


//sigma.setCapacity(params->nsnodes+2*WGDstree->nWGD);


*/


 //for (int i=params->nsnodes; i<WGDstree->nnodes; i++) {
 
for (int i=0; i<WGDstree->nWGD; i++) {
 
    //names[WGDstree->theWGD[i]->WGD_before->name]=names[WGDstree->theWGD[i]->WGD_after->name];

    //names[WGDstree->theWGD[i]->WGD_at->name]=names[WGDstree->theWGD[i]->WGD_after->name];

printf("\nvoila lele %d\n",i);

printf("\nvoila WGDstree->theWGD[i]->WGD_after->name %d\n",WGDstree->theWGD[i]->WGD_after->name);

printf("\net sa veleur %f\n",mu[WGDstree->theWGD[i]->WGD_after->name]);




mu[WGDstree->theWGD[i]->WGD_at->name]=mu[WGDstree->theWGD[i]->WGD_after->name];
mu[WGDstree->theWGD[i]->WGD_before->name]=mu[WGDstree->theWGD[i]->WGD_at->name];

printf("\net la vealeur at %f\n",mu[WGDstree->theWGD[i]->WGD_at->name]);
printf("\net la valeur before %f\n",mu[WGDstree->theWGD[i]->WGD_before->name]);



sigma[WGDstree->theWGD[i]->WGD_at->name]=sigma[WGDstree->theWGD[i]->WGD_after->name];
sigma[WGDstree->theWGD[i]->WGD_before->name]=sigma[WGDstree->theWGD[i]->WGD_at->name];


}




 // at the end 
 return new SpidirParams(names.size(), names, mu, sigma, alpha, beta);



}

/////////////////////////////////////////////////////////////////////////////////////
  
//try to write q function for the doomed 

void calcDoomTableWGD(Tree *WGDstree, float birthRate, float deathRate, 
                   double *doomtable)
{
    const double l = birthRate;
    const double u = deathRate;
    const double r = l - u;
    const double lu = l / u;
    double p0, p1;

    // get nodes in post order
    ExtendArray<Node*> nodes(0, WGDstree->nnodes);
    getTreePostOrder(WGDstree, &nodes);    

    indWGD=0;
    
    for (int i=0; i<WGDstree->nnodes; i++) {
        Node *node = nodes[i];
        
        if (node->isLeaf()) {
            doomtable[node->name] = -INFINITY;
        } else {          
            double prod = 0.0;            
            for (int j=0; j<node->nchildren; j++) {
                Node *child = node->children[j];


		//try to add something for WGD

		if (child->name == theWGD[indWGD]->WGD_at->name){

		 indWGD++;
		  
		 prod= (1-theWGD[indWGD]->lossProb) * doomtable[child->name]^2 + (theWGD[indWGD]->lossProb) * doomtable[child->name];


		}else{ 
                // compute u_t and P(t(c))
                const double t = child->dist;
                const double dc = exp(doomtable[child->name]);

                if (birthRate == deathRate) {
                    const double lt = l * t;
                    const double lt1 = 1.0 + lt;
                    p0 = lt / lt1;
                    p1 = 1.0 / lt1 / lt1;
                } else {
                    const double ert = exp(-r * t);
                    const double luert = l - u*ert;
                    p0 = (u - u * ert) / luert;
                    p1 = r*r * ert / luert / luert;
                }

                prod += log(p0 + dc * p1 / (1.0 - lu * p0 * dc));
		}
	    
	    }
            doomtable[node->name] = prod;
	}
    }
}


 //////////////////////////// end try



    

} // end of namespace spidir
