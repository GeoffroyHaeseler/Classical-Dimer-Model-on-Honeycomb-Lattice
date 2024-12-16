#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <cmath>
#include <string>
#include <time.h>

#include "measure.h"
#include "defs.h"
#include "next.h"
#include "nspFunctions.h"

extern int xB;
extern int yB;
extern char triB;
extern int nB;
extern int xA;
extern int yA;
extern char triA;
extern int nA;

extern int g_lattice[M][L][L][4];
extern double g_magneticField[2];

extern double g_mag_x;
extern double g_mag_z;

extern int g_Nspins;

extern int N_external_loop;
extern int wormSize;

extern double Etrotter;
extern double g_layerEnergy[M];
extern double g_diagonalEnergy[M];
extern double pA;

extern double OP_Moy;
extern double OP_SQ;
extern double OP_QA;

extern double g_energy;
extern double g_energy_sq;
extern double g_energy_qa;

extern double g_mag;
extern double g_magx;
extern double g_magy;
extern double g_mag_sq;
extern double g_mag_qa;

extern int g_length;

/// <summary>
/// Measure the lattice properties and store it in the evolution.txt file.
/// </summary>
/// <param name="T"></param> Lattice temperature
void latticeMeasure(double V, int Nconfig) {
	int ik;

	double ener = 0;
	double mag = 0;
	double magx = 0;
	double magy = 0;

	double test = 0;

	for (ik = 0; ik < M; ik++) {
		ener += calculateDiagonalEnergy(ik, V);
		ener += calculateLayerZeemanEnergy(ik);

		magx = calculateLayerMagX(ik);
		magy = calculateLayerMagY(ik);
		mag = calculateLayerMag(ik);
	}

	if (DoYouWantNeutronScatteringPlot) { updateNSP(); }
	if (DoYouWantTransitionGraphNSP) { transitionGraphNSP(); }
	if (DoYouWantFlavienPlot) { FlavienPlotUpdate(Nconfig); }

	if (DoYouWant2DPlot) { plot2DUpdate(); }

	g_energy += ener;
	g_energy_sq += ener * ener;
	g_energy_qa += ener * ener * ener * ener;
	g_magx += magx;
	g_magy += magy;
	g_mag += mag;
	g_mag_sq += mag * mag;
	g_mag_qa += mag * mag * mag * mag;

	double OP = calculateOrderParameter();
	OP_Moy += OP;
	OP_SQ += OP * OP;
	OP_QA += OP * OP * OP * OP;
}

/// <summary>
/// Average all the measures done on the lattice.
/// </summary>
void endMeasure(double T, double B, double V, double phi, int Nmeasurment) {

	g_energy /= Nmeasurment;
	g_energy /= g_Nspins;

	g_energy_sq /= Nmeasurment;
	g_energy_sq /= g_Nspins;
	g_energy_sq /= g_Nspins;

	g_energy_qa /= Nmeasurment;
	g_energy_qa /= g_Nspins;
	g_energy_qa /= g_Nspins;
	g_energy_qa /= g_Nspins;
	g_energy_qa /= g_Nspins;

	g_mag /= Nmeasurment;
	g_mag /= g_Nspins;

	g_mag_sq /= Nmeasurment;
	g_mag_sq /= g_Nspins;
	g_mag_sq /= g_Nspins;

	g_mag_qa /= Nmeasurment;
	g_mag_qa /= g_Nspins;
	g_mag_qa /= g_Nspins;
	g_mag_qa /= g_Nspins;
	g_mag_qa /= g_Nspins;

	g_magx /= Nmeasurment;
	g_magx /= g_Nspins;

	g_magy /= Nmeasurment;
	g_magy /= g_Nspins;

	OP_Moy /= Nmeasurment;

	OP_SQ /= Nmeasurment;

	OP_QA /= Nmeasurment;


	double Cb = g_energy_sq - g_energy * g_energy;
	Cb *= g_Nspins;
	Cb /= T * T;

	double EnergyBinderCumulent = 1.0 - g_energy_qa / (3.0 * g_energy_sq * g_energy_sq);

	double Chi = g_mag_sq - g_mag * g_mag;
	Chi *= g_Nspins;
	Chi /= T;

	double MagBinderCumulent = 1.0 - g_mag_qa / (3.0 * g_mag_sq * g_mag_sq);


	double newChi = OP_SQ - OP_Moy * OP_Moy;
	newChi *= g_Nspins;
	newChi /= T;
	
	double r3r3BinderCumulent = 1.0 - OP_QA / (3.0 * OP_SQ * OP_SQ);

	pA /= Nmeas;
	pA /= Neq;

	std::cout << std::endl;
	std::cout << "<E>  = " << g_energy << std::endl;
	std::cout << "<E2> = " << g_energy_sq << std::endl;
	std::cout << "<E4> = " << g_energy_qa << std::endl;
	std::cout << "Cb   = " << Cb << std::endl;
	std::cout << "BC   = " << EnergyBinderCumulent << std::endl;
	std::cout << std::endl;
	std::cout << "<m>  = " << g_mag << std::endl;
	std::cout << "<m2> = " << g_mag_sq << std::endl;
	std::cout << "Chi  = " << Chi << std::endl;
	std::cout << std::endl;
	std::cout << "OP_Moy = " << OP_Moy << std::endl;
	std::cout << "OP_Var = " << OP_SQ << std::endl;
	std::cout << "newChi = " << newChi << std::endl;
	//std::cout << std::endl;
	//std::cout << "pA  = " << pA << std::endl;
	std::cout << std::endl;



	//Save the equillibrium value in evolution.txt
	std::fstream myfile;
	myfile.open("evolution.txt", std::fstream::out | std::fstream::app);
	myfile << T << "\t" << phi << "\t" << B << "\t" << V << "\t" << Etrotter << "\t";
	myfile << "\t" << g_energy << "\t" << g_energy_sq << "\t" << Cb << "\t" << "\t";
	myfile << g_energy_qa << "\t" << EnergyBinderCumulent << "\t" << "\t";
	myfile << g_magx << "\t" << g_magy << "\t" << atan(g_magx / g_magy) << "\t";
	myfile << g_mag << "\t" << g_mag_sq << "\t" << Chi << "\t" << MagBinderCumulent << "\t\t";
	myfile << OP_Moy << "\t" << OP_SQ << "\t" << newChi << "\t" << r3r3BinderCumulent << std::endl;
	myfile.close();


	if (DoYouWantNeutronScatteringPlot) { endNSPMeasures(); }
	if (DoYouWantTransitionGraphNSP) { endNSPTransitionGraphMeasures(); }
	if (DoYouWantFlavienPlot) { endFlavienPlot(); }

	if (DoYouWant2DPlot) { end2DPlot(); }

	g_energy = 0;
	g_energy_sq = 0;
	g_energy_qa = 0;

	g_mag = 0;
	g_mag_sq = 0;
	g_mag_qa = 0;

	g_magx = 0;
	g_magy = 0;

	pA = 0;

	OP_Moy = 0;
	OP_SQ = 0;
	OP_QA = 0;
}



/// <summary>
/// Calculate the Zeeman energy of the kth layer.
/// </summary>
/// <param name="ix"></param>
/// <param name="iy"></param>
/// <param name="iz"></param>
/// <param name="itri"></param>
/// <param name="ispin"></param>
/// <returns></returns>
double calculateLayerEnergy(int k, double V) {
	double ener = 0;
	ener += calculateDiagonalEnergy(k, V);
	ener += calculateTrotterEnergy(k - 1);
	ener += calculateTrotterEnergy(k);
	ener += calculateLayerZeemanEnergy(k);
	return(ener);
}

/// <summary>
/// Calculate the energy du to the diagonal term of the hamiltonian (V term)
/// </summary>
/// <param name="k"></param>
/// <returns></returns>
double calculateDiagonalEnergy(int k, double V) {
	/*
	There is as many hexagones as up triangles.
	We loop over the up triangles (ix, iy, itri) and look for tha ssociated hexagone.

					/ \
				   /nx1\
				  /	ny1 \
				 / ntri1 \
	   _________/1_______2\_________
	   \	   /		   \       /
		\	  /			    \     /
		 \   /				 \   /
		  \ /				  \ /
		  /0\				  /0\
		 /   \				 /nx2\
		/ix,iy\				/ ny2 \
	   / itri  \		   / ntri2 \
	  /________2\_________/1________\
				 \		 /
				  \		/
				   \   /
					\ /

	*/

	int ix, iy, itri;
	double energy = 0;

	int nx1, nx2, ny1, ny2, ntri1, ntri2;
	int useless;

	for (ix = 0; ix < L; ix++) {
		for (iy = 0; iy < L; iy++) {
			for (itri = 0; itri < 4; itri++) {
				next(ix, iy, itri, 0, 1, 0, &nx1, &ny1, &ntri1, &useless);
				next(ix, iy, itri, 2, 1, 0, &nx2, &ny2, &ntri2, &useless);
				if (g_lattice[k][ix][iy][itri] == 0) {
					if (g_lattice[k][nx1][ny1][ntri1] == 2) {
						if (g_lattice[k][nx2][ny2][ntri2] == 1) {
							energy += V;
						}
					}
				}

				if (g_lattice[k][ix][iy][itri] == 2) {
					if (g_lattice[k][nx1][ny1][ntri1] == 1) {
						if (g_lattice[k][nx2][ny2][ntri2] == 0) {
							energy += V;
						}
					}
				}
			}
		}
	}
	return(energy);
}

/// <summary>
/// Calculate the energy du to ising coupling in trotter dimension between the kth and k+1th layer.
/// </summary>
/// <param name="ix"></param>
/// <param name="iy"></param>
/// <param name="iz"></param>
/// <param name="itri"></param>
/// <param name="ispin"></param>
/// <returns></returns>
double calculateTrotterEnergy(int k) {
	int ix, iy, itri;
	double ener = 0;

	int a = k + 1;
	if (a == M) { a = 0; }

	for (ix = 0; ix < L; ix++) {
		for (iy = 0; iy < L; iy++) {
			for (itri = 0; itri < 4; itri++) {
				if (g_lattice[k][ix][iy][itri] != g_lattice[a][ix][iy][itri]) { ener -= 1 * Etrotter; }
				else { ener += 3 * Etrotter; } // Etrotter < 0
			}
		}
	}
	return(ener);
}



/// <summary>
/// calculate the energy du to magnetic field interaction in all the kth layer
/// </summary>
/// <param name="k"></param>
/// <returns></returns>
double calculateLayerZeemanEnergy(int k) {
	int ix, iy, itri, ispin;
	double ener = 0;
	for (ix = 0; ix < L; ix++) {
		for (iy = 0; iy < L; iy++) {
			for (itri = 0; itri < 4; itri++) {
				for (ispin = 0; ispin < 3; ispin++) {
					ener += calculateMagEnergy(k, ix, iy, itri, ispin);
				}
			}
		}
	}
	return(ener);
}

/// <summary>
/// Calculate the Zeeman energy of the spin located at the given position.
/// </summary>
/// <param name="ix"></param>
/// <param name="iy"></param>
/// <param name="iz"></param>
/// <param name="itri"></param>
/// <param name="ispin"></param>
/// <returns></returns>
double calculateMagEnergy(int k, int ix, int iy, int itri, int ispin) {
	double ener = -(S[ispin][0] * g_magneticField[0] + S[ispin][1] * g_magneticField[1]);
	if (ispin == g_lattice[k][ix][iy][itri]) {//spin in
		ener = -ener;
	}
	return(ener / Sperp);
}



/// <summary>
/// Calculate the magnitisation along the X axis in the kth layer.
/// </summary>
/// <param name="k"></param>
/// <returns></returns>
double calculateLayerMagX(int k) {
	double magX = 0;
	int ix, iy, itri, ispin;
	for (ix = 0; ix < L; ix++) {
		for (iy = 0; iy < L; iy++) {
			for (itri = 0; itri < 4; itri++) {
				for (ispin = 0; ispin < 3; ispin++) {
					magX += calculateMagX(k, ix, iy, itri, ispin);
				}
			}
		}
	}
	return(magX);
}

/// <summary>
/// Calculate the magnetisation along the x axis of the spin located at the given position.
/// </summary>
/// <param name="ix"></param>
/// <param name="iy"></param>
/// <param name="iz"></param>
/// <param name="itri"></param>
/// <param name="ispin"></param>
/// <returns></returns>
double calculateMagX(int k, int ix, int iy, int itri, int ispin) {
	double mag_x = S[ispin][0]; //True if spin in
	if (ispin == g_lattice[k][ix][iy][itri]) {//spin out
		return (-mag_x);
	}
	return (mag_x);
}



/// <summary>
/// Calculate the magnitisation along the Y axis in the kth layer.
/// </summary>
/// <param name="k"></param>
/// <returns></returns>
double calculateLayerMagY(int k) {
	double magY = 0;
	int ix, iy, itri, ispin;
	for (ix = 0; ix < L; ix++) {
		for (iy = 0; iy < L; iy++) {
			for (itri = 0; itri < 4; itri++) {
				for (ispin = 0; ispin < 3; ispin++) {
					magY += calculateMagY(k, ix, iy, itri, ispin);
				}
			}
		}
	}
	return(magY);
}

/// <summary>
/// Calculate the magnetisation along the x axis of the spin located at the given position.
/// </summary>
/// <param name="ix"></param>
/// <param name="iy"></param>
/// <param name="iz"></param>
/// <param name="itri"></param>
/// <param name="ispin"></param>
/// <returns></returns>
double calculateMagY(int k, int ix, int iy, int itri, int ispin) {
	double mag_y = S[ispin][1]; //True if spin in
	if (ispin == g_lattice[k][ix][iy][itri]) {//spin out
		return (-mag_y);
	}
	return (mag_y);
}


/// <summary>
/// Calculate the topological sector in which is the system.
/// </summary>
/// <param name="k"></param>
/// <returns></returns>
int topologicalSector0(int k) {
	int ix;
	int sector0 = 0;
	for (ix = 0; ix < L; ix++) {
		if (g_lattice[k][ix][0][0] == 0) { sector0 += 1; }
		else { sector0 -= 1; }
		if (g_lattice[k][ix][0][1] == 0) { sector0 += 1; }
		else { sector0 -= 1; }
	}
	return(sector0);
}

/// <summary>
/// Calculate the topological sector in which is the system.
/// </summary>
/// <param name="k"></param>
/// <returns></returns>
int topologicalSector1(int k) {
	int ix;
	int sector1 = 0;
	for (ix = 0; ix < L; ix++) {
		if (g_lattice[k][ix][0][0] == 1) { sector1 += 1; }
		else { sector1 -= 1; }
		if (g_lattice[k][ix][0][1] == 1) { sector1 += 1; }
		else { sector1 -= 1; }
	}
	return(sector1);
}

/// <summary>
/// Calculate the topological sector in which is the system.
/// </summary>
/// <param name="k"></param>
/// <returns></returns>
int topologicalSector2(int k) {
	int ix;
	int sector2 = 0;
	for (ix = 0; ix < L; ix++) {
		if (g_lattice[k][ix][0][0] == 0) { sector2 += 1; }
		else { sector2 -= 1; }
		if (g_lattice[k][ix][0][1] == 0) { sector2 += 1; }
		else { sector2 -= 1; }
	}
	return(sector2);
}



/// <summary>
/// Add the current magnetisation to the file magnetizationHistogram.txt.
/// </summary>
void magnetizationHistogramMeasure() {
	double mag = 3.0 * calculateLayerMag(0) / (2.0 * g_Nspins * Sperp);
	std::fstream myfile;
	myfile.open("magnetizationHistogram.txt", std::fstream::out | std::fstream::app);
	myfile << mag << std::endl;
	myfile.close();
}

/// <summary>
/// Add the current r3r3 Order Paramter to the file Histogram.txt.
/// </summary>
void orderParameterHistogramMeasure() {
	double OP = calculateOrderParameter();
	double mag = 3.0 * calculateLayerMag(0) / (2.0 * g_Nspins * Sperp);
	std::fstream myfile;
	myfile.open("orderParameterHistogram.txt", std::fstream::out | std::fstream::app);
	myfile << OP << "\t" << mag << std::endl;
	myfile.close();
}



/// <summary>
/// Calculate the magnitisation along the B axis in the kth layer.
/// </summary>
/// <param name="k"></param>
/// <returns></returns>
double calculateLayerMag(int k) {
	double Bnorm = g_magneticField[0] * g_magneticField[0] + g_magneticField[1] * g_magneticField[1];
	if (Bnorm == 0) {
		double magX, magY;
		magX = calculateLayerMagX(k);
		magY = calculateLayerMagY(k);
		return(sqrt(magX * magX + magY * magY));
	}

	double mag = 0;
	int ix, iy, itri, ispin;
	for (ix = 0; ix < L; ix++) {
		for (iy = 0; iy < L; iy++) {
			for (itri = 0; itri < 4; itri++) {
				for (ispin = 0; ispin < 3; ispin++) {
					mag += calculateMag(k, ix, iy, itri, ispin);
				}
			}
		}
	}

	return(mag);
}

/// <summary>
/// Calculate the magnetisation along the magnetic field axis of the spin located at the given position.
/// </summary>
/// <param name="ix"></param>
/// <param name="iy"></param>
/// <param name="iz"></param>
/// <param name="itri"></param>
/// <param name="ispin"></param>
/// <returns></returns>
double calculateMag(int k, int ix, int iy, int itri, int ispin) {
	double Bnorm = g_magneticField[0] * g_magneticField[0] + g_magneticField[1] * g_magneticField[1];
	double mag = calculateMagX(k, ix, iy, itri, ispin) * g_magneticField[0] + calculateMagY(k, ix, iy, itri, ispin) * g_magneticField[1];

	return(mag / sqrt(Bnorm));
}



/// <summary>
/// Save the spin orientations in the lattice with the 'worm' in lattice.txt
/// </summary>
void save_lattice() {
	std::fstream myfile;
	std::remove("D:\\3A\\CF2\\lattice.txt");
	myfile.open("lattice.txt", std::fstream::out | std::fstream::app);

	int x = 0;
	int y = 0;
	int tri = 0;
	for (x = 0; x < L; x++) {
		for (y = 0; y < L; y++) {
			for (tri = 0; tri < 4; tri++) {

				myfile << g_lattice[0][x][y][tri] << std::endl;

			}
		}
	}
	myfile.close();

	//Save the initial conditions of the 'worm' in CI.txt
	myfile.open("CI.txt", std::fstream::out | std::fstream::app);

	myfile << xA << yA << triA << nA << std::endl;
	myfile << xB << yB << triB << nB << std::endl;
	myfile.close();

	std::cout << "Initial conditions are :" << std::endl;
	std::cout << "xA  = " << xA << std::endl;
	std::cout << "yA  = " << yA << std::endl;
	std::cout << "triA  = " << triA << std::endl;
	std::cout << "nA  = " << nA << std::endl;

	std::cout << "xB  = " << xB << std::endl;
	std::cout << "yB  = " << yB << std::endl;
	std::cout << "triB  = " << triB << std::endl;
	std::cout << "nB  = " << nB << std::endl;

	if (DoYouWantToSaveRef) { save_ref(); }

}

/// <summary>
/// Save the spin orientations in the lattice with the 'worm' as a reference
/// </summary>
void save_ref() {
	//Initialise the results files.
	std::remove("D:\\3A\\CF2\\lattice.txt");

	std::fstream myfile;
	myfile.open("D:\\3A\\CF2\\lattice.txt", std::fstream::out | std::fstream::app);

	int x = 0;
	int y = 0;
	int tri = 0;
	for (x = 0; x < L; x++) {
		for (y = 0; y < L; y++) {
			for (tri = 0; tri < 4; tri++) {

				myfile << g_lattice[0][x][y][tri] << std::endl;

			}
		}
	}
	myfile.close();

	//Save the initial conditions of the 'worm' in CI.txt
	myfile.open("D:\\3A\\CF2\\CI.txt", std::fstream::out | std::fstream::app);

	myfile << xA << yA << triA << nA << std::endl;
	myfile << xB << yB << triB << nB << std::endl;
	myfile.close();
}


/// <summary>
/// Return true if g_lattice[k] is in the star pahse topological sector. Return false otherwise.
/// </summary>
/// <param name="k"></param>
bool weAreInStarPhaseTopologicalSector(int k) {

	if (calculateLayerMagX(k) == 0 && calculateLayerMagY(k) == 0) { return(true); }
	return(false);
}