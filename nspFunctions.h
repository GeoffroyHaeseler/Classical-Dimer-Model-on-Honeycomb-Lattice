#include <complex>
using namespace std::complex_literals;

void nspInitialisation();
void updateNSP();
void endNSPMeasures();

void transitionGraphNSP();
void endNSPTransitionGraphMeasures();

void FlavienInitialisation();
void FlavienPlotUpdate(int Nconfig);
void endFlavienPlot();

void plot2DInitialisation();
void plot2DUpdate();
void end2DPlot();

double calculateOrderParameter();

double dotProduct(double A[3], const double B[3]);
double dotProduct(double A[3], double B[3]);
std::complex<double> dotProduct(std::complex<double> A[3], std::complex<double> B[3]);
double norm(double A[3]);
std::complex<double> norm(std::complex<double> A[3]);