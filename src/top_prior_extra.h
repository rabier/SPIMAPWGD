#ifndef SPIDIR_TOP_PRIOR_EXTRA_H
#define SPIDIR_TOP_PRIOR_EXTRA_H

#include "Tree.h"


namespace spidir {

extern "C" {


int inumHistories(int ngenes);

double numHistories(int ngenes);

int inumTopologyHistories(Tree *tree);

double numTopologyHistories(Tree *tree);


double birthDeathTreePrior2(Tree *tree, Tree *stree, int *recon, 
                          int *events, float birth, float death,
			    double *doomtable, float q);


}

} // namespace spidir

#endif // SPIDIR_TOP_PRIOR_EXTRA_H
