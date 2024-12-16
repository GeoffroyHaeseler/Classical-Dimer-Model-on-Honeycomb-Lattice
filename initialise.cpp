#include <iostream>
#include <fstream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <cmath>

#include "initialise.h"
#include "addAWorm.h"
#include "measure.h"
#include "nspFunctions.h"
#include "defs.h"

using namespace std::complex_literals;

extern int g_lattice[M][L][L][4];
extern bool g_TypeAHexagone[L][L][4];
extern bool g_TypeCHexagone[L][L][4];
extern double z0;
extern double z1;
extern double z2;
extern double kT;
extern double proba[3][3][3][3][3][3][3];
extern double g_magneticField[2];

extern double Etrotter;
extern double pA;
extern double OP_Moy;
extern double OP_SQ;
extern double OP_QA;

extern double g_energy;
extern double g_energy_sq;
extern double g_energy_qa;
extern double g_mag;
extern double g_mag_sq;
extern double g_mag_qa;

extern std::complex<double> g_phaseOP[L][L][4][3];
extern double g_SpinPosition[L][L][4][4][3];

extern int g_Nspins;

/// <summary>
/// Initialise the lattice and all the relevant arrays and values before the code begin
/// </summary>
void initialise() {
	//Initialise the random seed.
	srand(time(NULL));
	double A = (double)rand() / RAND_MAX;

	//Initialise the statistical averages
	g_energy = 0;
	g_energy_sq = 0;
	g_energy_qa = 0;
	g_mag = 0;
	g_mag_sq = 0;
	g_mag_qa = 0;

	g_Nspins = 12 * L * L;

	OP_Moy = 0;
	OP_SQ = 0;
	OP_QA = 0;
	pA = 0;

	//Initialise the lattice.
	lattice2Initialisation();
	intialiseZ3SymmetryBreaking();

	spinPositionIntialisation();

	if (DoYouWantNeutronScatteringPlot || DoYouWantTransitionGraphNSP) { nspInitialisation(); }
	if (DoYouWantFlavienPlot) { FlavienInitialisation(); }
	if (DoYouWant2DPlot) { plot2DInitialisation(); }

	//Initialise the results files.
	std::remove("cube.ply");
	std::remove("results.ply");
	std::remove("evolution.txt");
	std::remove("magnetizationHistogram.txt");
	std::remove("orderParameterHistogram.txt");

	std::fstream myfile;
	myfile.open("evolution.txt", std::fstream::out | std::fstream::app);
	myfile << "T (J)\tPhi (rad)\tB0 (J)\tV0 (J)\tEtrotter\t\t<E>\t<E²>\tCb\t\t<E4>\tBinder\t\t<mx>\t<my>\tangle\t<m>\t<m²>\tChi\tMagBinderCumulent\t\tOP_Moy\tOP_Var\tnewChi\tr3r3BinderCumulent\t\tL = " << L << "\tNspins = " << g_Nspins << "\tNmeas = " << Nmeas << "\tNeq = " << Neq << "\tTi = " << Ti << "\tTf = " << Tf << "\tdT = " << dT << std::endl;
	myfile.close();
}

/// <summary>
/// initialise the magnetic field (in the ux, uy, uz base, cf next.cpp).
/// </summary>
/// <param name="theta"></param>
/// <param name="phi"></param>
void magFieldInitialisation(double B, double phi) {
	g_magneticField[0] = B * sin(phi);
	g_magneticField[1] = B * cos(phi);
	//std::cout << g_magneticField[0] << std::endl << g_magneticField[1] << std::endl;
}

/// <summary>
/// Intialise the lattice. We begin in the columnar phase.
/// </summary>
void latticeInitialisation() {
	std::cout << "The lattice is initialised with a columnar phase\n";
	int ik, ix, iy, itri;
	for (ik = 0; ik < M; ik++) {
		for (ix = 0; ix < L; ix++) {
			for (iy = 0; iy < L; iy++) {
				for (itri = 0; itri < 4; itri++) {
					g_lattice[ik][ix][iy][itri] = 0;
				}
			}
		}
	}
}

/// <summary>
/// Intialise the lattice. We begin in the star phase.
/// </summary>
void lattice2Initialisation() {
	std::cout << "The lattice is initialised with a star phase (2)\n";

	int ik, ix, iy;
	for (ik = 0; ik < M; ik++) {
		for (ix = 0; ix < L / 3.0; ix++) {
			for (iy = 0; iy < L / 3.0; iy++) {
				g_lattice[ik][3 * ix + 0][3 * iy + 0][0] = 0;
				g_lattice[ik][3 * ix + 0][3 * iy + 0][1] = 1;
				g_lattice[ik][3 * ix + 0][3 * iy + 0][2] = 0;
				g_lattice[ik][3 * ix + 0][3 * iy + 0][3] = 2;

				g_lattice[ik][3 * ix + 1][3 * iy + 0][0] = 2;
				g_lattice[ik][3 * ix + 1][3 * iy + 0][1] = 0;
				g_lattice[ik][3 * ix + 1][3 * iy + 0][2] = 2;
				g_lattice[ik][3 * ix + 1][3 * iy + 0][3] = 1;

				g_lattice[ik][3 * ix + 2][3 * iy + 0][0] = 1;
				g_lattice[ik][3 * ix + 2][3 * iy + 0][1] = 2;
				g_lattice[ik][3 * ix + 2][3 * iy + 0][2] = 1;
				g_lattice[ik][3 * ix + 2][3 * iy + 0][3] = 0;



				g_lattice[ik][3 * ix + 0][3 * iy + 1][0] = 0;
				g_lattice[ik][3 * ix + 0][3 * iy + 1][1] = 1;
				g_lattice[ik][3 * ix + 0][3 * iy + 1][2] = 0;
				g_lattice[ik][3 * ix + 0][3 * iy + 1][3] = 2;

				g_lattice[ik][3 * ix + 1][3 * iy + 1][0] = 2;
				g_lattice[ik][3 * ix + 1][3 * iy + 1][1] = 0;
				g_lattice[ik][3 * ix + 1][3 * iy + 1][2] = 2;
				g_lattice[ik][3 * ix + 1][3 * iy + 1][3] = 1;

				g_lattice[ik][3 * ix + 2][3 * iy + 1][0] = 1;
				g_lattice[ik][3 * ix + 2][3 * iy + 1][1] = 2;
				g_lattice[ik][3 * ix + 2][3 * iy + 1][2] = 1;
				g_lattice[ik][3 * ix + 2][3 * iy + 1][3] = 0;



				g_lattice[ik][3 * ix + 0][3 * iy + 2][0] = 0;
				g_lattice[ik][3 * ix + 0][3 * iy + 2][1] = 1;
				g_lattice[ik][3 * ix + 0][3 * iy + 2][2] = 0;
				g_lattice[ik][3 * ix + 0][3 * iy + 2][3] = 2;

				g_lattice[ik][3 * ix + 1][3 * iy + 2][0] = 2;
				g_lattice[ik][3 * ix + 1][3 * iy + 2][1] = 0;
				g_lattice[ik][3 * ix + 1][3 * iy + 2][2] = 2;
				g_lattice[ik][3 * ix + 1][3 * iy + 2][3] = 1;

				g_lattice[ik][3 * ix + 2][3 * iy + 2][0] = 1;
				g_lattice[ik][3 * ix + 2][3 * iy + 2][1] = 2;
				g_lattice[ik][3 * ix + 2][3 * iy + 2][2] = 1;
				g_lattice[ik][3 * ix + 2][3 * iy + 2][3] = 0;
			}
		}
	}
}

/// <summary>
/// Intialise the lattice. We begin in the star phase.
/// </summary>
void lattice3Initialisation() {
	std::cout << "The lattice is initialised with a star phase (3)\n";
	int ik, ix, iy;
	for (ik = 0; ik < M; ik++) {
		for (ix = 0; ix < L / 3.0; ix++) {
			for (iy = 0; iy < L / 3.0; iy++) {
				g_lattice[ik][3 * ix + 0][3 * iy + 0][0] = 2;
				g_lattice[ik][3 * ix + 0][3 * iy + 0][1] = 0;
				g_lattice[ik][3 * ix + 0][3 * iy + 0][2] = 2;
				g_lattice[ik][3 * ix + 0][3 * iy + 0][3] = 1;

				g_lattice[ik][3 * ix + 1][3 * iy + 0][0] = 1;
				g_lattice[ik][3 * ix + 1][3 * iy + 0][1] = 2;
				g_lattice[ik][3 * ix + 1][3 * iy + 0][2] = 1;
				g_lattice[ik][3 * ix + 1][3 * iy + 0][3] = 0;

				g_lattice[ik][3 * ix + 2][3 * iy + 0][0] = 0;
				g_lattice[ik][3 * ix + 2][3 * iy + 0][1] = 1;
				g_lattice[ik][3 * ix + 2][3 * iy + 0][2] = 0;
				g_lattice[ik][3 * ix + 2][3 * iy + 0][3] = 2;



				g_lattice[ik][3 * ix + 0][3 * iy + 1][0] = 2;
				g_lattice[ik][3 * ix + 0][3 * iy + 1][1] = 0;
				g_lattice[ik][3 * ix + 0][3 * iy + 1][2] = 2;
				g_lattice[ik][3 * ix + 0][3 * iy + 1][3] = 1;

				g_lattice[ik][3 * ix + 1][3 * iy + 1][0] = 1;
				g_lattice[ik][3 * ix + 1][3 * iy + 1][1] = 2;
				g_lattice[ik][3 * ix + 1][3 * iy + 1][2] = 1;
				g_lattice[ik][3 * ix + 1][3 * iy + 1][3] = 0;

				g_lattice[ik][3 * ix + 2][3 * iy + 1][0] = 0;
				g_lattice[ik][3 * ix + 2][3 * iy + 1][1] = 1;
				g_lattice[ik][3 * ix + 2][3 * iy + 1][2] = 0;
				g_lattice[ik][3 * ix + 2][3 * iy + 1][3] = 2;



				g_lattice[ik][3 * ix + 0][3 * iy + 2][0] = 2;
				g_lattice[ik][3 * ix + 0][3 * iy + 2][1] = 0;
				g_lattice[ik][3 * ix + 0][3 * iy + 2][2] = 2;
				g_lattice[ik][3 * ix + 0][3 * iy + 2][3] = 1;

				g_lattice[ik][3 * ix + 1][3 * iy + 2][0] = 1;
				g_lattice[ik][3 * ix + 1][3 * iy + 2][1] = 2;
				g_lattice[ik][3 * ix + 1][3 * iy + 2][2] = 1;
				g_lattice[ik][3 * ix + 1][3 * iy + 2][3] = 0;

				g_lattice[ik][3 * ix + 2][3 * iy + 2][0] = 0;
				g_lattice[ik][3 * ix + 2][3 * iy + 2][1] = 1;
				g_lattice[ik][3 * ix + 2][3 * iy + 2][2] = 0;
				g_lattice[ik][3 * ix + 2][3 * iy + 2][3] = 2;
			}
		}
	}
}

/// <summary>
/// Intialise the lattice. We begin in the star phase.
/// </summary>
void lattice4Initialisation() {
	std::cout << "The lattice is initialised with a star phase (4)\n";
	int ik, ix, iy;
	for (ik = 0; ik < M; ik++) {
		for (ix = 0; ix < L / 3.0; ix++) {
			for (iy = 0; iy < L / 3.0; iy++) {
				g_lattice[ik][3 * ix + 0][3 * iy + 0][0] = 1;
				g_lattice[ik][3 * ix + 0][3 * iy + 0][1] = 2;
				g_lattice[ik][3 * ix + 0][3 * iy + 0][2] = 1;
				g_lattice[ik][3 * ix + 0][3 * iy + 0][3] = 0;

				g_lattice[ik][3 * ix + 1][3 * iy + 0][0] = 0;
				g_lattice[ik][3 * ix + 1][3 * iy + 0][1] = 1;
				g_lattice[ik][3 * ix + 1][3 * iy + 0][2] = 0;
				g_lattice[ik][3 * ix + 1][3 * iy + 0][3] = 2;

				g_lattice[ik][3 * ix + 2][3 * iy + 0][0] = 2;
				g_lattice[ik][3 * ix + 2][3 * iy + 0][1] = 0;
				g_lattice[ik][3 * ix + 2][3 * iy + 0][2] = 2;
				g_lattice[ik][3 * ix + 2][3 * iy + 0][3] = 1;



				g_lattice[ik][3 * ix + 0][3 * iy + 1][0] = 1;
				g_lattice[ik][3 * ix + 0][3 * iy + 1][1] = 2;
				g_lattice[ik][3 * ix + 0][3 * iy + 1][2] = 1;
				g_lattice[ik][3 * ix + 0][3 * iy + 1][3] = 0;

				g_lattice[ik][3 * ix + 1][3 * iy + 1][0] = 0;
				g_lattice[ik][3 * ix + 1][3 * iy + 1][1] = 1;
				g_lattice[ik][3 * ix + 1][3 * iy + 1][2] = 0;
				g_lattice[ik][3 * ix + 1][3 * iy + 1][3] = 2;

				g_lattice[ik][3 * ix + 2][3 * iy + 1][0] = 2;
				g_lattice[ik][3 * ix + 2][3 * iy + 1][1] = 0;
				g_lattice[ik][3 * ix + 2][3 * iy + 1][2] = 2;
				g_lattice[ik][3 * ix + 2][3 * iy + 1][3] = 1;



				g_lattice[ik][3 * ix + 0][3 * iy + 2][0] = 1;
				g_lattice[ik][3 * ix + 0][3 * iy + 2][1] = 2;
				g_lattice[ik][3 * ix + 0][3 * iy + 2][2] = 1;
				g_lattice[ik][3 * ix + 0][3 * iy + 2][3] = 0;

				g_lattice[ik][3 * ix + 1][3 * iy + 2][0] = 0;
				g_lattice[ik][3 * ix + 1][3 * iy + 2][1] = 1;
				g_lattice[ik][3 * ix + 1][3 * iy + 2][2] = 0;
				g_lattice[ik][3 * ix + 1][3 * iy + 2][3] = 2;

				g_lattice[ik][3 * ix + 2][3 * iy + 2][0] = 2;
				g_lattice[ik][3 * ix + 2][3 * iy + 2][1] = 0;
				g_lattice[ik][3 * ix + 2][3 * iy + 2][2] = 2;
				g_lattice[ik][3 * ix + 2][3 * iy + 2][3] = 1;
			}
		}
	}
}



/// <summary>
/// initialise g_TypeAHexagone[L][L][4];
/// </summary>
void intialiseZ3SymmetryBreaking() {
	for (int ix = 0; ix < L / 3.0; ix++) {
		for (int iy = 0; iy < L; iy++) {
			g_TypeAHexagone[3 * ix + 0][iy][0] = true;
			g_TypeAHexagone[3 * ix + 0][iy][1] = false;
			g_TypeAHexagone[3 * ix + 0][iy][2] = true;
			g_TypeAHexagone[3 * ix + 0][iy][3] = false;

			g_TypeAHexagone[3 * ix + 1][iy][0] = false;
			g_TypeAHexagone[3 * ix + 1][iy][1] = true;
			g_TypeAHexagone[3 * ix + 1][iy][2] = false;
			g_TypeAHexagone[3 * ix + 1][iy][3] = false;

			g_TypeAHexagone[3 * ix + 2][iy][0] = false;
			g_TypeAHexagone[3 * ix + 2][iy][1] = false;
			g_TypeAHexagone[3 * ix + 2][iy][2] = false;
			g_TypeAHexagone[3 * ix + 2][iy][3] = true;
		}
	}

	for (int ix = 0; ix < L / 3.0; ix++) {
		for (int iy = 0; iy < L; iy++) {
			g_TypeCHexagone[3 * ix + 0][iy][0] = false;
			g_TypeCHexagone[3 * ix + 0][iy][1] = false;
			g_TypeCHexagone[3 * ix + 0][iy][2] = false;
			g_TypeCHexagone[3 * ix + 0][iy][3] = true;

			g_TypeCHexagone[3 * ix + 1][iy][0] = true;
			g_TypeCHexagone[3 * ix + 1][iy][1] = false;
			g_TypeCHexagone[3 * ix + 1][iy][2] = true;
			g_TypeCHexagone[3 * ix + 1][iy][3] = false;

			g_TypeCHexagone[3 * ix + 2][iy][0] = false;
			g_TypeCHexagone[3 * ix + 2][iy][1] = true;
			g_TypeCHexagone[3 * ix + 2][iy][2] = false;
			g_TypeCHexagone[3 * ix + 2][iy][3] = false;
		}
	}

}

/// <summary>
/// Calculate the position of each spin in the lattice and their associate phase
/// </summary>
void spinPositionIntialisation() {
	//Initialise the Real Lattice
	int iX, iY, itri, ispin;
	for (iX = 0; iX < L; iX++) {//We go through the whole real lattice
		for (iY = 0; iY < L; iY++) {

			for (itri = 0; itri < 4; itri++) {
				for (ispin = 0; ispin < 4; ispin++) {
					g_SpinPosition[iX][iY][itri][ispin][0] = 2.0 * iX * ReaL[0][0] + iY * (2.0 * ReaL[1][0] - ReaL[0][0]);
					g_SpinPosition[iX][iY][itri][ispin][1] = 2.0 * iX * ReaL[0][1] + iY * (2.0 * ReaL[1][1] - ReaL[0][1]);
					g_SpinPosition[iX][iY][itri][ispin][2] = 2.0 * iX * ReaL[0][2] + iY * (2.0 * ReaL[1][2] - ReaL[0][2]);
				}
			}

			for (ispin = 0; ispin < 4; ispin++) {
				g_SpinPosition[iX][iY][0][ispin][0] += 0;
				g_SpinPosition[iX][iY][0][ispin][1] += 0;
				g_SpinPosition[iX][iY][0][ispin][2] += 0;

				g_SpinPosition[iX][iY][1][ispin][0] += ReaL[0][0];
				g_SpinPosition[iX][iY][1][ispin][1] += ReaL[0][1];
				g_SpinPosition[iX][iY][1][ispin][2] += ReaL[0][2];

				g_SpinPosition[iX][iY][2][ispin][0] += ReaL[0][0] + ReaL[1][0];
				g_SpinPosition[iX][iY][2][ispin][1] += ReaL[0][1] + ReaL[1][1];
				g_SpinPosition[iX][iY][2][ispin][2] += ReaL[0][2] + ReaL[1][2];

				g_SpinPosition[iX][iY][3][ispin][0] += ReaL[1][0];
				g_SpinPosition[iX][iY][3][ispin][1] += ReaL[1][1];
				g_SpinPosition[iX][iY][3][ispin][2] += ReaL[1][2];
			}

			for (itri = 0; itri < 4; itri++) {
				g_SpinPosition[iX][iY][itri][0][0] += 0.5 * ReaL[1][0];
				g_SpinPosition[iX][iY][itri][0][1] += 0.5 * ReaL[1][1];
				g_SpinPosition[iX][iY][itri][0][2] += 0.5 * ReaL[1][2];

				g_SpinPosition[iX][iY][itri][1][0] += 0;
				g_SpinPosition[iX][iY][itri][1][1] += 0;
				g_SpinPosition[iX][iY][itri][1][2] += 0;

				g_SpinPosition[iX][iY][itri][2][0] += 0.5 * ReaL[0][0];
				g_SpinPosition[iX][iY][itri][2][1] += 0.5 * ReaL[0][1];
				g_SpinPosition[iX][iY][itri][2][2] += 0.5 * ReaL[0][2];

				g_SpinPosition[iX][iY][itri][3][0] += 0.5 * ReaL[2][0];
				g_SpinPosition[iX][iY][itri][3][1] += 0.5 * ReaL[2][1];
				g_SpinPosition[iX][iY][itri][3][2] += 0.5 * ReaL[2][2];
			}

		}
	}

	const std::complex<double> i(0, 1);
	double Q[3] = { 2 * pi2 / (3.0 * NT), 0, 0 };

	for (iX = 0; iX < L; iX++) {//We go through the whole real lattice
		for (iY = 0; iY < L; iY++) {
			for (itri = 0; itri < 4; itri++) {
				for (ispin = 0; ispin < 3; ispin++) {
					g_phaseOP[iX][iY][itri][ispin] = std::exp(i * dotProduct(Q, g_SpinPosition[iX][iY][itri][ispin]));

				}
			}
		}
	}

}

/// <summary>
/// Initialise the probability matrix. Return true if the T>Tc
/// </summary>
void probaInitialisation(double T, double B, double V, double phi) {
	//initialise the magnetic field.
	magFieldInitialisation(B, phi);

	//Initialise the different probability involved
	double mu0, mu1, mu2;
	mu0 = +4.0 * sqrt(2) * B *  cos(phi) / 3.0;
	mu1 = -2.0 * sqrt(2) * B * (cos(phi) + sqrt(3) * sin(phi)) / 3.0;
	mu2 = -2.0 * sqrt(2) * B * (cos(phi) - sqrt(3) * sin(phi)) / 3.0;

	kT = T;

	//Boltzman weights for the case of purely Zeeman interaction
	z0 = exp(mu0 / T);
	z1 = exp(mu1 / T);
	z2 = exp(mu2 / T);

	double Kt = 0; //Ising constant between layers
	//double Kt = -0.5 * log(tanh(g / (kT * M))); //Ising constant between layers
	Etrotter = -T * Kt; //>0 if spins are different, the energy is increased.

	double plaquetteWeight[3];

	plaquetteWeight[0] = exp(-0.0 * V / T); //exp(-bV)
	plaquetteWeight[1] = exp(-1.0 * V / T);
	plaquetteWeight[2] = exp(-2.0 * V / T);

	double deltaWeight[3];
	deltaWeight[0] = exp(-1.0 * delta / T);
	deltaWeight[1] = exp(-0.0 * delta / T);
	deltaWeight[2] = exp(+1.0 * delta / T);

	double w0, w1, w2;

	int Vcount0, Vcount1, Vcount2, deltacount0, deltacount1, deltacount2;
	int Ktcount1, Ktcount2;
	for (Vcount0 = 0; Vcount0 < 3; Vcount0++) {
		for (Vcount1 = 0; Vcount1 < 3; Vcount1++) {
			for (Vcount2 = 0; Vcount2 < 3; Vcount2++) {
				for (deltacount0 = 0; deltacount0 < 3; deltacount0++) {
					for (deltacount1 = 0; deltacount1 < 3; deltacount1++) {
						for (deltacount2 = 0; deltacount2 < 3; deltacount2++) {
							//update the Boltzman weights with the diagonal energy (V term)
							w0 = z0 * plaquetteWeight[Vcount0] * deltaWeight[deltacount0];
							w1 = z1 * plaquetteWeight[Vcount1] * deltaWeight[deltacount1];
							w2 = z2 * plaquetteWeight[Vcount2] * deltaWeight[deltacount2];

							proba[Vcount0][Vcount1][Vcount2][deltacount0][deltacount1][deltacount2][0] = w0 / (w0 + w1 + w2);	//P(0)
							proba[Vcount0][Vcount1][Vcount2][deltacount0][deltacount1][deltacount2][1] = w1 / (w0 + w1 + w2);	//P(1)
							proba[Vcount0][Vcount1][Vcount2][deltacount0][deltacount1][deltacount2][2] = w2 / (w0 + w1 + w2);	//P(2)
							//std::cout << proba[Vcount0][Vcount1][Vcount2][deltacount0][deltacount1][deltacount2][0] << "\t" << proba[Vcount0][Vcount1][Vcount2][deltacount0][deltacount1][deltacount2][1] << "\t" << proba[Vcount0][Vcount1][Vcount2][deltacount0][deltacount1][deltacount2][2] << "\n";
						}
					}
				}
			}
		}
	}
	std::cout << "Lattice at T = " << T << "J, B = " << B << "J and V = " << V << "J is now being measured..." << std::endl;
	std::cout << "In this case the quantum coupling constant is : Kt = " << Kt << std::endl;
}
