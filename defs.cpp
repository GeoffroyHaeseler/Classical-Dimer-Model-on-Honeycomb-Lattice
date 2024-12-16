#include "defs.h"
#include <complex>

using namespace std::complex_literals;

int g_lattice[M][L][L][4];
int g_latticeTransitionGraphLattice[L][L][4][3];
bool g_TypeAHexagone[L][L][4];
bool g_TypeCHexagone[L][L][4];
double z0;
double z1;
double z2;
double kT;
double g_magneticField[2];

double OP_Moy;
double OP_SQ;
double OP_QA;

int g_Nspins;

//Coupling constants
double Etrotter;//Interaction energy coupling layers.
double proba[3][3][3][3][3][3][3];


double pA;

//thermodynamical properties
double g_magx;
double g_magy;
double g_mag;
double g_mag_sq;
double g_mag_qa;

double g_energy;
double g_energy_sq;
double g_energy_qa;

double g_diagonalEnergy[M];
double g_layerEnergy[M];

//addAWorm intial defect neighbours
int xA;
int yA;
char triA;
int nA;

int xB;
int yB;
char triB;
int nB;

int xC;
int yC;
char triC;
int nC;

//NSP values

std::complex<double> g_phaseOP[L][L][4][3];
double g_SpinPosition[L][L][4][4][3];
double g_Q[4 * xdiv][4 * ydiv][3];
double g_QHat[4 * xdiv][4 * ydiv][3];
double g_AllInAllOutOrder[4 * xdiv][4 * ydiv];
double g_SecondeFragmentation[4 * xdiv][4 * ydiv];
double g_UnpolarisedNSP[4 * xdiv][4 * ydiv];
double g_PseudoSpinNSP[4 * xdiv][4 * ydiv];
double g_InPlaneNSP[4 * xdiv][4 * ydiv];
double g_BraggsPeak[4 * xdiv][4 * ydiv];
double g_delta[3];
double g_qlim;
int g_nspCount;

//Flavien's plot
double gF_Q3D[R][angle][angle][3];
double gF_QHat3D[R][angle][angle][3];
double gF_Qperp3D[R][angle][angle][3];
double gF_QperpHat3D[R][angle][angle][3];
double gF_Q2D[R][angle][3];
double gF_QHat2D[R][angle][3];
double gF_2DPlot[R];
double gF_plot[R][2];
double gF_magnetic;
double gF_nuclear;
double gF_QlimDw;
double gF_QlimUp;
int gF_layer;