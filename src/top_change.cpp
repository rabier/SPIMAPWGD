/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gene tree search functions

   modified for WGD by Charles-Elie Rabier
    2011-2013
=============================================================================*/

#include "common.h"
#include "Tree.h"


namespace spidir {


//=============================================================================
// Nearest Neighbor Interchange Topology Proposal

/*

    Proposes a new tree using Nearest Neighbor Interchange
       
       Branch for NNI is specified by giving its two incident nodes (node1 and 
       node2).  Change specifies which  subtree of node1 will be swapped with
       the uncle.  See figure below.

         node2
        /     \
      nodeb    node1
               /  \
         nodea     * 

*/
void performNni(Tree *tree, Node *nodea, Node *nodeb)
{
    Node *node1 = nodea->parent;
    Node *node2 = nodeb->parent;
    
    // assert that node1 and node2 are incident to the same branch
    assert(node1->parent == node2 ||
           node2->parent == node1);
    
    // find child indexes
    int a = (node1->children[0] == nodea) ? 0 : 1;
    assert(node1->children[a] == nodea);

    int b = (node2->children[0] == nodeb) ? 0 : 1;
    assert(node2->children[b] == nodeb);
    
    // swap parent pointers
    nodea->parent = node2;
    nodeb->parent = node1;
    
    // swap child pointers
    node2->children[b] = nodea;
    node1->children[a] = nodeb;
}






void proposeRandomNni(Tree *tree, Node **a, Node **b)
{
    // find edges for NNI
    int choice;
    do {
        choice = irand(tree->nnodes);
    } while (tree->nodes[choice]->isLeaf() || 
             tree->nodes[choice]->parent == NULL);
    
    Node *node1 = tree->nodes[choice];
    Node *node2 = tree->nodes[choice]->parent;
    *a = node1->children[irand(2)];
    *b = (node2->children[0] == node1) ? node2->children[1] :
                                         node2->children[0];
    assert((*a)->parent->parent == (*b)->parent);
}



//========================================== SubtreeSlide is the LOCAL CHANGE of LARGET and SIMON

/*      noded
          \
         node2
        /     \
      nodec    node1
               /  \
         nodea     nodeb 

*/




  void performSubtreeSlide(Tree *tree, float *mratio, float *mstarratio)
{

  // find an  edge such as 
  //the nodes must not be leaves
  //but a node of the edge can be the root
  int choice;
 
  do {
    choice = irand(tree->nnodes);
  } while (tree->nodes[choice]->isLeaf() || 
	   tree->nodes[choice]->parent == NULL) ;

  
  Node *node1 = tree->nodes[choice];
  Node *node2 = tree->nodes[choice]->parent;    
  Node  *nodec =(node2->children[0] == node1) ? node2->children[1] :
    node2->children[0];
  Node *noded=node2->parent;
  Node *nodea = node1->children[irand(2)];
  //chose randomly a children of node 1
  //we call it nodea
  //nodea will be a an extremity  for the local change algorithm

  //nodeb is the other children
  Node  *nodeb =(node1->children[0] == nodea) ? node1->children[1] :
                                 node1->children[0];
  float u1=frand();
  float lambda=0.2;
  float u2=frand();
  float x,y,xstar,ystar;

  if (frand()<0.5){

    //nodec is the other extremity for the local change        
    float m = nodea->dist + node1->dist + nodec->dist;

    //m is the distance nodea nodec
    *mratio=m;
    float mstar=m*exp(lambda*(u1-0.5));
    *mstarratio=mstar;
    x=nodea->dist;
    y=nodea->dist+node1->dist;    
    
    if (frand()<0.5){
      // node2  does not move
      //ie node1 moves
      xstar=u2*mstar;
      ystar=y*mstar/m;
      
      if (xstar<ystar){
	//topology does not change
	nodea->dist=xstar;
	node1->dist=ystar-xstar;
	nodec->dist=mstar-ystar;
      }else{
	//topology changes
	if (node1->children[0]==nodeb){
	  node1->children[1]=nodec;
	}else{
	  node1->children[0]=nodec;
	}
	
	if (node2->children[0]==node1){
	  node2->children[1]=nodea;
	}else{
	  node2->children[0]=nodea;
	}

	nodec->parent=node1; 
	nodea->parent=node2;	
	nodea->dist=ystar;
	node1->dist=xstar-ystar;
	nodec->dist=mstar-xstar;
      }

    }else{
      //node1 does not move
      //node2 moves
      ystar=u2*mstar;
      xstar=x*mstar/m;
	
      if (xstar<ystar){
	//topology does not change
	nodea->dist=xstar;
	node1->dist=ystar-xstar;
	nodec->dist=mstar-ystar;
      }else{
	//topology changes
	
	if (node1->children[0]==nodeb){
	  node1->children[1]=nodec;
	}else{
	  node1->children[0]=nodec;
	}
	
	nodec->parent=node1;
	
	if (node2->children[0]==node1){
	  node2->children[1]=nodea;
	}else{
	  node2->children[0]=nodea;
	}

	nodea->parent=node2;	
	nodea->dist=ystar;
	node1->dist=xstar-ystar;
	nodec->dist=mstar-xstar;
	
      }
    }
    
  }else{

    //noded is the other extremity for the local change
    float m = nodea->dist + node1->dist + node2->dist;
    //m is the distance nodea noded
    *mratio=m;
    float mstar=m*exp(lambda*(u1-0.5));
    *mstarratio=mstar;
    x=nodea->dist;
    y=nodea->dist+node1->dist;
    
    if (frand()<0.5){
      // node2  does not move
      //ie node 1 moves      
      xstar=u2*mstar;
      ystar=y*mstar/m;
      
      if (xstar<ystar){
	//topology does not change
	nodea->dist=xstar;
	node1->dist=ystar-xstar;
	node2->dist=mstar-ystar;
      }else{

	if (node2->parent==NULL){
	  //topology can not change
	  //just change the length 
	  
	  do{
	    mstar=m*exp(lambda*(frand()-0.5));
	    xstar=frand()*mstar;
	    ystar=mstar;} while (xstar>ystar);
	  
	  nodea->dist=xstar;
	  node1->dist=ystar-xstar;
	  
	}else{
	  //topology changes	
	  
	  if (node2->children[0]==nodec){
	    node2->children[1]=nodea;
	  }else{
	    node2->children[0]=nodea;
	  } 

	  nodea->parent=node2;

	  if (node1->children[0]==nodeb){
	    node1->children[1]=node2;
	  }else{
	    node1->children[0]=node2;
	  } 
	   
	  node2->parent=node1;
	  if (noded->children[0]==node2){
	    noded->children[0]=node1;
	  }else{
	    noded->children[1]=node1;
	  } 

	  node1->parent=noded;
	  nodea->dist=ystar;
	  node2->dist=xstar-ystar;
	  node1->dist=mstar-xstar;
	}
      }


    }else{
      //node1 does not move
      //ie node 2 moves
      //node2 has to have a parent otherwise it is not possible
      //so, in case node2 does not have any parent , we perform subtreeslide  again 
      if (node2->parent==NULL){
	performSubtreeSlide(tree,mratio,mstarratio);
      }else{
	
	//printf("\n nod2parent\n");
	ystar=u2*mstar;
	xstar=x*mstar/m;
	
	if (xstar<ystar){
	  //topology does not change
	  nodea->dist=xstar;
	  node1->dist=ystar-xstar;
	  node2->dist=mstar-ystar;

	}else{
	  //topology changes
	 
	  if (node2->children[0]==nodec){
	    node2->children[1]=nodea;
	  }else{
	    node2->children[0]=nodea;}
	  
	  if (node1->children[0]==nodeb){
	    node1->children[1]=node2;
	  }else{
	    node1->children[0]=node2;}
	  
	  if (noded->children[0]==node2){
	    noded->children[0]=node1;
	  }else{
	    noded->children[1]=node1;
	  }

	  nodea->parent=node2;
	  node2->parent=node1;
	  node1->parent=noded;
	  nodea->dist=ystar;
	  node2->dist=xstar-ystar;
	  node1->dist=mstar-xstar;

	}
      }
    }
  }
  
}


  //======================================================================
  //Change just a little bite the length of one edge of the tree

 void performBranchLength(Tree *tree, float *mratio, float *mstarratio)
{
  // find a node which is not the leaf
  //we will change the length of the edge above this node

  int choice;
  do {
    choice = irand(tree->nnodes);   
  } while (tree->nodes[choice]->parent == NULL) ;
       
  Node *node1 = tree->nodes[choice];
  float lambda=0.2;
  float u1=frand();
  float m = node1->dist ;
  float mstar=m*exp(lambda*(u1-0.5));
  node1->dist = mstar;

  *mratio=m;
  *mstarratio=mstar;


}


//=============================================================================
// Subtree pruning and regrafting (SPR)


/*
    a = subtree
    e = newpos
    
    BEFORE
            ....
        f         d
       /           \
      c             e
     / \           ...
    a   b
   ... ...

    AFTER

        f         d
       /           \
      b             c
     ...           / \
                  a   e
                 ... ...

    Requirements:
    1. a (subtree) is not root or children of root
    2. e (newpos) is not root, a, descendant of a, c (parent of a), or 
       b (sibling of a)
    3. tree is binary

*/
void performSpr(Tree *tree, Node *subtree, Node *newpos)
{
    Node *a = subtree;
    Node *e = newpos;

    Node *c = a->parent;
    Node *f = c->parent;
    const int bi = (c->children[0] == a) ? 1 : 0;
    Node *b = c->children[bi];
    const int ci = (f->children[0] == c) ? 0 : 1;
    Node *d = e->parent;
    const int ei = (d->children[0] == e) ? 0 : 1;

    d->children[ei] = c;
    c->children[bi] = e;
    f->children[ci] = b;
    b->parent = f;
    c->parent = d;
    e->parent = c;
}

/*
    What if e == f  (also equivalent to NNI) this is OK

    BEFORE
    
          d
         / \
        e  ...
       / \
      c  ...         
     / \           
    a   b
   ... ...

    AFTER
          d
         / \
        c
       / \
      a   e
     ... / \
        b  ...
       ...
       
  What if d == f  (also equivalent to NNI) this is OK
  
    BEFORE
          
        f
       / \
      c   e
     / \  ...
    a   b
   ... ...

    AFTER
          
        f
       / \
      b   c  
     ... / \ 
        a   e
       ... ...  
*/



/*
    Requirements:
    1. a (subtree) is not root or children of root
    2. e (newpos) is not root, a, descendant of a, c (parent of a), or 
       b (sibling of a)
    3. tree is binary
*/
void proposeRandomSpr(Tree *tree, Node **subtree, Node **newpos)
{
    assert(tree->nnodes >= 5);

    // find subtree (a) to cut off (any node that is not root or child of root)
    int choice;
        do {
        choice = irand(tree->nnodes);
    } while (tree->nodes[choice]->parent == NULL ||
             tree->nodes[choice]->parent->parent == NULL);
    Node *a = tree->nodes[choice];
    *subtree = a;
    
    // find sibling (b) of a
    const Node *c = a->parent;
    const int bi = (c->children[0] == a) ? 1 : 0;
    const Node *b = c->children[bi];
    
    // choose newpos (e)
    Node *e = NULL;
    do {
        choice = irand(tree->nnodes);
        e = tree->nodes[choice];
        
        // test if e is a valid choice
        if (e->parent == NULL || e == a || e == c || e == b)
            continue;
        
        // also test if e is a descendent of a
        bool under_a = false;
        for (Node *ptr = e->parent; ptr != NULL; ptr = ptr->parent) {
            if (ptr == a) {
                under_a = true;
                break;
            }
        }            
        
        if (under_a)
            continue;
        
        break;
    } while (true);
    *newpos = e;
}


bool validSpr(Tree *tree, const Node *subtree, const Node *newpos)
{
    const Node *a = subtree;
    const Node *e = newpos;
    
    // find sibling (b) of a
    const Node *c = a->parent;
    const int bi = (c->children[0] == a) ? 1 : 0;
    const Node *b = c->children[bi];

    // test if a is a valid choice
    if (a->parent == NULL || a->parent->parent == NULL)
        return false;

    // test if e is a valid choice
    if (e->parent == NULL || e == a || e == c || e == b) {
        //printf("ERROR: a=%d, e=%d, b=%d, c=%d\n", a->name, e->name, 
        //       b->name, c->name);
        return false;
    }
        
    // also test if e is a descendent of a
    for (Node *ptr = e->parent; ptr != NULL; ptr = ptr->parent)
        if (ptr == a)
            return false;
        
    return true;
}

} // namespace spidir
