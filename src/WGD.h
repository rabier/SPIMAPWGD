#define SPIDIR_NEWICK_H

#include <stdlib.h>
#include <string>

#include "Tree.h"
#include "phylogeny.h"
#include "model_params.h"

namespace spidir {

SpeciesTree  *removeWGDnodes(SpeciesTree *tree);
void extendRateParamToWGDnodes(SpidirParams *params, SpeciesTree *WGDstree);

}
