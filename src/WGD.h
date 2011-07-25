#define SPIDIR_NEWICK_H

#include <stdlib.h>
#include <string>

#include "Tree.h"
#include "phylogeny.h"
#include "model_params.h"

namespace spidir {

//void removeWGDnodes(Tree *tree);
SpeciesTree  *removeWGDnodes(SpeciesTree *tree);
 SpidirParams *paramforWGD(SpidirParams *params, SpeciesTree *tree);

}
