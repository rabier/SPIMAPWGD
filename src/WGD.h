#define SPIDIR_NEWICK_H

#include <stdlib.h>
#include <string>

#include "Tree.h"
#include "phylogeny.h"
#include "model_params.h"

namespace spidir {

SpeciesTree  *removeWGDnodes(SpeciesTree *tree);
void extendRateParamToWGDnodes(SpidirParams *params, SpeciesTree *WGDstree);
int  WGDreconcile_onebranch(SpeciesTree *WGDstree,  int *recon_noWGD, int *events, Node *nodeWGDafter,Node *node,int *best);
void WGDreconReset(SpeciesTree *WGDstree,  int *recon, int *recon_noWGD, int *events, Node *nodeWGDafter,Node *node,int *best);

}
