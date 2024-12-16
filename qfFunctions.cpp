#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#include "qfFunctions.h"
#include "defs.h"
#include "addAWorm.h"
#include "measure.h"
#include "loop.h"

extern int g_lattice[M][L][L][4];

extern double g_layerEnergy[M];
extern double g_diagonalEnergy[M];

extern double pA;

bool addAQuantumWorm() {

	addAWorm(0); // add a Loop in g_lattice

	return(true);
}
