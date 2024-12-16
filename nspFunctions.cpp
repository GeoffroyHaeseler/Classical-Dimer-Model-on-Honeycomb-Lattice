#include <iostream>
#include <fstream>
#include <stdlib.h>     /* srand, rand */
#include <cmath>
#include <math.h>
#include <complex.h>
#include <string>
#include <time.h>

#include "nspFunctions.h"
#include "defs.h"

extern int g_lattice[M][L][L][4];
extern int g_latticeTransitionGraphLattice[L][L][4][3];
extern int g_Nspins;

extern double g_SpinPosition[L][L][4][4][3];
extern std::complex<double> g_phaseOP[L][L][4][3];
extern double g_Q[4 * xdiv][4 * ydiv][3];
extern double g_QHat[4 * xdiv][4 * ydiv][3];
extern double g_AllInAllOutOrder[4 * xdiv][4 * ydiv];
extern double g_SecondeFragmentation[4 * xdiv][4 * ydiv];
extern double g_UnpolarisedNSP[4 * xdiv][4 * ydiv];
extern double g_PseudoSpinNSP[4 * xdiv][4 * ydiv];
extern double g_InPlaneNSP[4 * xdiv][4 * ydiv];
extern double g_BraggsPeak[4 * xdiv][4 * ydiv];
extern double g_delta[3];

extern double g_qlim;
extern int g_nspCount;

extern double gF_Q3D[R][angle][angle][3];
extern double gF_QHat3D[R][angle][angle][3];
extern double gF_Qperp3D[R][angle][angle][3];
extern double gF_QperpHat3D[R][angle][angle][3];
extern double gF_Q2D[R][angle][3];
extern double gF_QHat2D[R][angle][3];
extern double gF_2DPlot[R];
extern double gF_plot[R][2];
extern double gF_magnetic;
extern double gF_nuclear;
extern double gF_QlimDw;
extern double gF_QlimUp;
extern int gF_layer;

using namespace std::complex_literals;

/// <summary>
/// Initialise g_Q and g_QHat
/// </summary>
void nspInitialisation() {
	g_Nspins = 12 * L * L;//(3*(2L)² )
	g_nspCount = 0;

	//Initialise the Reciprocal Lattice
	int q0, q1;
	double Q0, Q1;
	double normsq;
	for (q0 = 0; q0 < 4 * xdiv; q0++) {//We go through the whole reciprocal lattice
		for (q1 = 0; q1 < 4 * ydiv; q1++) {
            Q0 = (double(q0) - 2.0 * xdiv) * xmult / (double(xdiv));
			Q1 = (double(q1) - 2.0 * ydiv) * ymult / (double(ydiv));

			g_Q[q0][q1][0] = Q0 * RecL[0][0] + Q1 * RecL[1][0];
			g_Q[q0][q1][1] = Q0 * RecL[0][1] + Q1 * RecL[1][1];
			g_Q[q0][q1][2] = Q0 * RecL[0][2] + Q1 * RecL[1][2];

			normsq = norm(g_Q[q0][q1]);

			if (normsq < 0.00001) {
				g_QHat[q0][q1][0] = 0;
				g_QHat[q0][q1][1] = 0;
			}
			else {
				g_QHat[q0][q1][0] = g_Q[q0][q1][0] / (sqrt(normsq));
				g_QHat[q0][q1][1] = g_Q[q0][q1][1] / (sqrt(normsq));
			}
			g_QHat[q0][q1][2] = 0;

			g_AllInAllOutOrder[q0][q1] = 0;
			g_SecondeFragmentation[q0][q1] = 0;
			g_UnpolarisedNSP[q0][q1] = 0;
			g_PseudoSpinNSP[q0][q1] = 0;
			g_InPlaneNSP[q0][q1] = 0;
			g_BraggsPeak[q0][q1] = 0;

		}
	}

	g_qlim = std::abs(g_Q[0][ydiv][0]);
}

/// <summary>
/// Update the values of the NSP (bragg's peaks, unpolarised and polarised neutron scattering).
/// </summary>
void updateNSP() {
	g_nspCount++;
	std::cout << "Enter update nsp function n " << g_nspCount << ".";
	time_t time_stamp_1 = time(NULL);

	int q0, q1;
	int iX, iY, itri, ispin;
	double sign;
	double ScalarProduct;
	const std::complex<double> i(0, 1);
	std::complex<double> phase = 0;
	std::complex<double> spinProjectedOnQ[3] = { 0, 0, 0 };
	std::complex<double> braggspeak = 0;
	std::complex<double> unpolarised[3] = { 0, 0, 0 };
	std::complex<double> SecondeFragmentation[2] = { 0, 0 };
	std::complex<double> AllInAllOutOrder[2] = { 0, 0 };
	double normsq;
	for (q0 = 0; q0 < 4 * xdiv; q0++) {//loop over the reciprocal lattice
		for (q1 = 0; q1 < 4 * ydiv; q1++) {
			normsq = norm(g_Q[q0][q1]);
			if (sqrt(normsq) < g_qlim) {
				for (iX = 0; iX < L; iX++) {//loop over the real lattice
					for (iY = 0; iY < L; iY++) {
						for (itri = 0; itri < 4; itri++) {

							/*//Switch from a square lattice (cf next.cpp) to a trigonal lattice (cf.    )
							X = 2 * i0;
							Y = 2 * i1;
							if (itri == 1 || itri == 2) { X++; }
							if (itri == 2 || itri == 3) {
								Y++;
								X -= i1;
								if (X < 0) { X += 2 * L; }
							}
							if (itri == 0 || itri == 1) {
								X -= i1;
								if (X < 0) { X += 2 * L; }
							}*/


							for (ispin = 0; ispin < 3; ispin++) {
								if (g_lattice[0][iX][iY][itri] == ispin) { sign = -1; }
								else { sign = +1; }

								braggspeak += std::exp(i  * dotProduct(g_Q[q0][q1], g_SpinPosition[iX][iY][itri][ispin]));

								phase = std::exp(i * sign * dotProduct(g_Q[q0][q1], g_SpinPosition[iX][iY][itri][ispin]));

								ScalarProduct = dotProduct(g_QHat[q0][q1], S[ispin]);
								spinProjectedOnQ[0] = ScalarProduct * g_QHat[q0][q1][0];
								spinProjectedOnQ[1] = ScalarProduct * g_QHat[q0][q1][1];
								spinProjectedOnQ[2] = ScalarProduct * g_QHat[q0][q1][2];

								//the whole spin is taken into account (unpolarised neutrons)
								unpolarised[0] += sign * (S[ispin][0] - spinProjectedOnQ[0]) * phase;
								unpolarised[1] += sign * (S[ispin][1] - spinProjectedOnQ[1]) * phase;
								unpolarised[2] += sign * (S[ispin][2] - spinProjectedOnQ[2]) * phase;
								SecondeFragmentation[0] += (sign - 1.0 / 3.0) * (S[ispin][0] - spinProjectedOnQ[0]) * phase;
								SecondeFragmentation[1] += (sign - 1.0 / 3.0) * (S[ispin][1] - spinProjectedOnQ[1]) * phase;
								AllInAllOutOrder[0] += (1.0 / 3.0) * (S[ispin][0] - spinProjectedOnQ[0]) * phase;
								AllInAllOutOrder[1] += (1.0 / 3.0) * (S[ispin][1] - spinProjectedOnQ[1]) * phase;
							}
						}
					}
				}//end of the loop on the real lattice
			}
			g_BraggsPeak[q0][q1] = std::abs(std::conj(braggspeak) * braggspeak);
			braggspeak = 0;

			//The whole spin is considered (Unpolarised Neutron Scattering)
			g_UnpolarisedNSP[q0][q1] += std::abs(norm(unpolarised));

			//Only the perpendicular component of the spin is taken into account (Polarised Neutron Scattering)
			g_PseudoSpinNSP[q0][q1] += std::abs(std::conj(unpolarised[2]) * unpolarised[2]);

			//Only the in kagome plane component of the spin is taken into account (Polarised Neutron Scattering)
			g_InPlaneNSP[q0][q1] += std::abs(std::conj(unpolarised[0]) * unpolarised[0] + std::conj(unpolarised[1]) * unpolarised[1]);

			//Only the dimer fragment of the kagome plane component of the spin is taken into account (Polarised Neutron Scattering)
			g_SecondeFragmentation[q0][q1] += std::abs(std::conj(SecondeFragmentation[0]) * SecondeFragmentation[0] + std::conj(SecondeFragmentation[1]) * SecondeFragmentation[1]);

			//Only the dimer fragment of the kagome plane component of the spin is taken into account (Polarised Neutron Scattering)
			g_AllInAllOutOrder[q0][q1] += std::abs(std::conj(AllInAllOutOrder[0]) * AllInAllOutOrder[0] + std::conj(AllInAllOutOrder[1]) * AllInAllOutOrder[1]);

			unpolarised[0] = 0;
			unpolarised[1] = 0;
			unpolarised[2] = 0;
			SecondeFragmentation[0] = 0;
			SecondeFragmentation[1] = 0;
			AllInAllOutOrder[0] = 0;
			AllInAllOutOrder[1] = 0;
		}
	}
	std::cout << " It took " << time(NULL) - time_stamp_1 << " seconds to finish the NSP update" << std::endl;
}

/// <summary>
/// Edit the NSP.txt file which contains the data requiered to plot the figure with the plot_AHC matlab's function.
/// </summary>
void endNSPMeasures() {

	std::cout << g_nspCount << " configurations have been used to calculated the NSP." << std::endl << std::endl;

	std::remove("NSP.txt");

	std::fstream myfile;
	myfile.open("NSP.txt", std::fstream::out | std::fstream::app);

	int q0, q1;
	double normsq;
	double maxUnpolarisedNSP = 0;
	double maxPseudoSpinNSP = 0;
	double maxInPlaneNSP = 0;
	double maxBraggsPeak = 0;

	//go through the whole reciprocal lattice to print the result in the NSP.txt file
	for (q0 = 0; q0 < 4 * xdiv; q0++) {
		for (q1 = 0; q1 < 4 * xdiv; q1++) {

			g_Q[q0][q1][0] *= 0.666666666;
			g_Q[q0][q1][1] *= 0.666666666;
			g_Q[q0][q1][2] *= 0.666666666;

			normsq = norm(g_Q[q0][q1]);

			if (sqrt(normsq) > g_qlim) {
				myfile << g_Q[q0][q1][0] << " " << g_Q[q0][q1][1] << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
			}

			else {
				g_AllInAllOutOrder[q0][q1] /= double(g_Nspins) * g_nspCount;
				g_SecondeFragmentation[q0][q1] /= double(g_Nspins) * g_nspCount;
				g_UnpolarisedNSP[q0][q1] /= double(g_Nspins) * g_nspCount;
				g_PseudoSpinNSP[q0][q1]  /= double(g_Nspins) * g_nspCount;
				g_InPlaneNSP[q0][q1]	 /= double(g_Nspins) * g_nspCount;
				g_BraggsPeak[q0][q1]	 /= double(g_Nspins) * g_nspCount;
				myfile << g_Q[q0][q1][0] << " " << g_Q[q0][q1][1] << " " << g_UnpolarisedNSP[q0][q1] << " ";
				myfile << g_PseudoSpinNSP[q0][q1] << " " << g_InPlaneNSP[q0][q1] << " " << g_BraggsPeak[q0][q1] << " " << g_SecondeFragmentation[q0][q1] << " " << g_AllInAllOutOrder[q0][q1] << std::endl;
			}
			g_AllInAllOutOrder[q0][q1] = 0;
			g_SecondeFragmentation[q0][q1] = 0;
			g_UnpolarisedNSP[q0][q1] = 0;
			g_PseudoSpinNSP[q0][q1] = 0;
			g_InPlaneNSP[q0][q1] = 0;
			g_BraggsPeak[q0][q1] = 0;
		}
	}

	myfile.close();
	g_nspCount = 0;
}



/// <summary>
/// Update the values of the NSP for transition graphs (bragg's peaks, unpolarised and polarised neutron scattering).
/// </summary>
void transitionGraphNSP() {
	g_nspCount++;
	std::cout << "Enter update nsp function n " << g_nspCount << ".";
	time_t time_stamp_1 = time(NULL);

	int q0, q1;
	int iX, iY, itri, ispin;
	double sign;
	double ScalarProduct;
	const std::complex<double> i(0, 1);
	std::complex<double> phase = 0;
	std::complex<double> spinProjectedOnQ[3] = { 0, 0, 0 };
	std::complex<double> braggspeak = 0;
	std::complex<double> unpolarised[3] = { 0, 0, 0 };
	double normsq;
	for (q0 = 0; q0 < 4 * xdiv; q0++) {//loop over the reciprocal lattice
		for (q1 = 0; q1 < 4 * ydiv; q1++) {
			normsq = norm(g_Q[q0][q1]);
			if (sqrt(normsq) < g_qlim) {
				for (iX = 0; iX < L; iX++) {//loop over the real lattice
					for (iY = 0; iY < L; iY++) {
						for (itri = 0; itri < 4; itri++) {


							for (ispin = 0; ispin < 3; ispin++) {
								sign = 0.5 * g_latticeTransitionGraphLattice[iX][iY][itri][ispin];
								//std::cout << sign << "\t";
								braggspeak += std::exp(i * dotProduct(g_Q[q0][q1], g_SpinPosition[iX][iY][itri][ispin]));

								phase = std::exp(i * sign * dotProduct(g_Q[q0][q1], g_SpinPosition[iX][iY][itri][ispin]));

								ScalarProduct = dotProduct(g_QHat[q0][q1], S[ispin]);
								spinProjectedOnQ[0] = ScalarProduct * g_QHat[q0][q1][0];
								spinProjectedOnQ[1] = ScalarProduct * g_QHat[q0][q1][1];
								spinProjectedOnQ[2] = ScalarProduct * g_QHat[q0][q1][2];

								//the whole spin is taken into account (unpolarised neutrons)
								unpolarised[0] += sign * (S[ispin][0] - spinProjectedOnQ[0]) * phase;
								unpolarised[1] += sign * (S[ispin][1] - spinProjectedOnQ[1]) * phase;
								unpolarised[2] += sign * (S[ispin][2] - spinProjectedOnQ[2]) * phase;
							}
							//std::cout << "\n";
						}
					}
				}//end of the loop on the real lattice
			}
			g_BraggsPeak[q0][q1] = std::abs(std::conj(braggspeak) * braggspeak);
			braggspeak = 0;

			//The whole spin is considered (Unpolarised Neutron Scattering)
			g_UnpolarisedNSP[q0][q1] += std::abs(norm(unpolarised));


			//Only the perpendicular component of the spin is taken into account (Polarised Neutron Scattering)
			g_PseudoSpinNSP[q0][q1] += std::abs(std::conj(unpolarised[2]) * unpolarised[2]);


			//Only the in kagome plane component of the spin is taken into account (Polarised Neutron Scattering)
			g_InPlaneNSP[q0][q1] += std::abs(std::conj(unpolarised[0]) * unpolarised[0] + std::conj(unpolarised[1]) * unpolarised[1]);

			unpolarised[0] = 0;
			unpolarised[1] = 0;
			unpolarised[2] = 0;
		}
	}
	std::cout << " It took " << time(NULL) - time_stamp_1 << " seconds to finish the NSP update" << std::endl;
}

/// <summary>
/// Edit the NSP.txt file which contains the data requiered to plot the NSP. Transition Graphs version.
/// </summary>
void endNSPTransitionGraphMeasures() {

	std::cout << g_nspCount << " configurations have been used to calculated the NSP." << std::endl << std::endl;

	std::remove("NSP.txt");

	std::fstream myfile;
	myfile.open("NSP.txt", std::fstream::out | std::fstream::app);

	int q0, q1;
	double normsq;
	double maxUnpolarisedNSP = 0;
	double maxPseudoSpinNSP  = 0;
	double maxInPlaneNSP     = 0;
	double maxBraggsPeak     = 0;

	//go through the whole reciprocal lattice to print the result in the NSP.txt file
	for (q0 = 0; q0 < 4 * xdiv; q0++) {
		for (q1 = 0; q1 < 4 * xdiv; q1++) {

			g_Q[q0][q1][0] *= 0.666666666;
			g_Q[q0][q1][1] *= 0.666666666;
			g_Q[q0][q1][2] *= 0.666666666;

			normsq = norm(g_Q[q0][q1]);

			if (sqrt(normsq) > g_qlim) {
				myfile << g_Q[q0][q1][0] << " " << g_Q[q0][q1][1] << " " << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
			}

			else {
				g_UnpolarisedNSP[q0][q1] /= double(g_Nspins) * Nconfiguration;
				g_PseudoSpinNSP[q0][q1]  /= double(g_Nspins) * Nconfiguration;
				g_InPlaneNSP[q0][q1]     /= double(g_Nspins) * Nconfiguration;
				g_BraggsPeak[q0][q1]     /= double(g_Nspins) * Nconfiguration;
				myfile << g_Q[q0][q1][0] << " " << g_Q[q0][q1][1] << " " << g_UnpolarisedNSP[q0][q1] << " ";
				myfile << g_PseudoSpinNSP[q0][q1] << " " << g_InPlaneNSP[q0][q1] << " " << g_BraggsPeak[q0][q1] << std::endl;
			}

			g_UnpolarisedNSP[q0][q1] = 0;
			g_PseudoSpinNSP[q0][q1]  = 0;
			g_InPlaneNSP[q0][q1]     = 0;
			g_BraggsPeak[q0][q1]     = 0;
		}
	}

	myfile.close();
	g_nspCount = 0;
}



/// <summary>
/// Create the spheres in reciprocal space on which FlavienPlotUpdate will integrated on.
/// </summary>
void FlavienInitialisation() {
	//gF_QlimDw = 1.0 * pi2 / NT;;
	//gF_QlimUp = 2.0 * pi2 / NT;

	gF_QlimDw = 0 * pi2 / NT;;
	gF_QlimUp = 3 * pi2 / NT;

	gF_layer = 0;

	int i0, i1;
	int ir, iphi, itheta;
	double r;
	double phi, theta;

	//intialise 2D parameters
	for (ir = 0; ir < R; ir++) {
		for (iphi = 0; iphi < angle; iphi++) {
			phi = pi2 * iphi / angle;
			r = gF_QlimDw + (gF_QlimUp - gF_QlimDw) * ir / R;

			gF_Q2D[ir][iphi][0] = r * cos(phi);
			gF_Q2D[ir][iphi][1] = r * sin(phi);
			gF_Q2D[ir][iphi][2] = 0;

			if (r < 0.00001) {
				gF_QHat2D[ir][iphi][0] = 0;
				gF_QHat2D[ir][iphi][1] = 0;
				gF_QHat2D[ir][iphi][2] = 0;
			}
			else {
				gF_QHat2D[ir][iphi][0] = cos(phi);
				gF_QHat2D[ir][iphi][1] = sin(phi);
				gF_QHat2D[ir][iphi][2] = 0;
			}

		}
	}

	//intialise 3D parameters
	for (ir = 0; ir < R; ir++) {
		for (itheta = 0; itheta < angle; itheta++) {
			for (iphi = 0; iphi < angle; iphi++) {
				r = gF_QlimDw + (gF_QlimUp - gF_QlimDw) * ir / double(R);
				theta = acos(1.0 - 2.0 * itheta / double(angle));
				phi = pi2 * iphi / double(angle);

				gF_Q3D[ir][iphi][itheta][0] = r * sin(theta) * cos(phi);
				gF_Q3D[ir][iphi][itheta][1] = r * sin(theta) * sin(phi);
				gF_Q3D[ir][iphi][itheta][2] = r * cos(theta);

				gF_Qperp3D[ir][iphi][itheta][0] = r * sin(theta) * cos(phi);
				gF_Qperp3D[ir][iphi][itheta][1] = r * sin(theta) * sin(phi);
				gF_Qperp3D[ir][iphi][itheta][2] = 0;

				if (r == 0) {
					gF_QHat3D[ir][iphi][itheta][0] = 0;
					gF_QHat3D[ir][iphi][itheta][1] = 0;
					gF_QHat3D[ir][iphi][itheta][2] = 0;
					gF_QperpHat3D[ir][iphi][itheta][0] = 0;
					gF_QperpHat3D[ir][iphi][itheta][1] = 0;
					gF_QperpHat3D[ir][iphi][itheta][2] = 0;
				}
				else {
					gF_QHat3D[ir][iphi][itheta][0] = sin(theta) * cos(phi);
					gF_QHat3D[ir][iphi][itheta][1] = sin(theta) * sin(phi);
					gF_QHat3D[ir][iphi][itheta][2] = cos(theta);

					gF_QperpHat3D[ir][iphi][itheta][0] = sin(theta) * cos(phi);
					gF_QperpHat3D[ir][iphi][itheta][1] = sin(theta) * sin(phi);
					gF_QperpHat3D[ir][iphi][itheta][2] = 0;
				}
			}
		}
	}

	g_delta[0] = ReaL[2][0];
	g_delta[1] = ReaL[2][1];
	g_delta[2] = ReaL[2][2];
}

/// <summary>
/// Integrate the intensity of the neutron scattering on the sphere of |Q| = cste.
/// </summary>
void FlavienPlotUpdate(int Nconfig) {
	time_t initial_time = time(NULL);
	std::cout << "Flavien's plot update n" << gF_layer << " is being calculated...";

	int ir, j;
	gF_nuclear = 0;
	gF_magnetic = 0;

	int iX, iY, itri, ispin;

	double r;
	double sintheta, costheta;
	double phi, sinphi, cosphi;
	double rng;

	double Q3D[3] = { 0, 0, 0 };
	double Qperp3D[3] = { 0, 0, 0 };

	double QHat3D[3] = { 0, 0, 0 };
	double QperpHat3D[3] = { 0, 0, 0 };

	double sign;
	double ScalarProduct;
	double spinProjectedOnQ[3] = { 0, 0, 0 };
	double delta[3];
	const std::complex<double> i(0, 1);
	std::complex<double> phase = 0;
	std::complex<double> unpolarised[3] = { 0, 0, 0 };

	std::complex<double> nuclear = 0;


	for (ir = 0; ir < R; ir++) {

		if (integratedOnASphere) {
			int itheta, iphi;
			for (itheta = 0; itheta < angle; itheta++) {
				for (iphi = 0; iphi < angle; iphi++) {

					r = gF_QlimDw + (gF_QlimUp - gF_QlimDw) * ir / double(R);

					costheta = 1.0 - 2.0 * itheta / double(angle);
					sintheta = 2.0 * sqrt(itheta * (angle - itheta)) / double(angle);;

					sinphi = sin(pi2 * iphi / double(angle));
					cosphi = cos(pi2 * iphi / double(angle));

					Q3D[0] = r * sintheta * cosphi;
					Q3D[1] = r * sintheta * sinphi;
					Q3D[2] = r * costheta;

					Qperp3D[0] = r * sintheta * cosphi;
					Qperp3D[1] = r * sintheta * sinphi;
					Qperp3D[2] = 0;


					if (r == 0) {
						QHat3D[0] = 0;
						QHat3D[1] = 0;
						QHat3D[2] = 0;
						QperpHat3D[0] = 0;
						QperpHat3D[1] = 0;
						QperpHat3D[2] = 0;
					}
					else {
						QHat3D[0] = sintheta * cosphi;
						QHat3D[1] = sintheta * sinphi;
						QHat3D[2] = costheta;

						QperpHat3D[0] = sintheta * cosphi;
						QperpHat3D[1] = sintheta * sinphi;
						QperpHat3D[2] = 0;
					}

					for (iX = 0; iX < L; iX++) {//loop over the real lattice
						for (iY = 0; iY < L; iY++) {
							for (itri = 0; itri < 4; itri++) {

								for (ispin = 0; ispin < 4; ispin++) {


									if (ferromagneticPeak) {
										if (ispin == 3) { sign = 1.0; } //the non kagome spins are out.
										else { sign = 1.0 / 9.0; }
										//sign = -1 if the spin is out, else it is +1. (because the direction of the spin in S is in,
										//and then must be flipped in function of the direction of the actual spin)

										//phase du to the spin position
										phase = std::exp(i * dotProduct(Q3D, g_SpinPosition[iX][iY][itri][ispin]));

										ScalarProduct = dotProduct(QHat3D, S[ispin]);
										spinProjectedOnQ[0] = ScalarProduct * QHat3D[0];
										spinProjectedOnQ[1] = ScalarProduct * QHat3D[1];
										spinProjectedOnQ[2] = ScalarProduct * QHat3D[2];

										//the whole spin is taken into account (unpolarised neutrons)
										unpolarised[0] += sign * phase;
										unpolarised[1] += sign * phase;
										unpolarised[2] += sign * phase;

										nuclear += sign * phase;
									}

									else {
										if (ispin == 3) { sign = -1; } //the non kagome spins are out.
										else if (ispin == g_lattice[0][iX][iY][itri]) { sign = -1; }
										else { sign = +1; }
										//sign = -1 if the spin is out, else it is +1. (because the direction of the spin in S is in,
										//and then must be flipped in function of the direction of the actual spin)

										//phase du to the spin position
										phase = std::exp(i * dotProduct(Q3D, g_SpinPosition[iX][iY][itri][ispin]));

										ScalarProduct = dotProduct(QHat3D, S[ispin]);
										spinProjectedOnQ[0] = ScalarProduct * QHat3D[0];
										spinProjectedOnQ[1] = ScalarProduct * QHat3D[1];
										spinProjectedOnQ[2] = ScalarProduct * QHat3D[2];

										//the whole spin is taken into account (unpolarised neutrons)
										unpolarised[0] += sign * (S[ispin][0] - spinProjectedOnQ[0]) * phase;
										unpolarised[1] += sign * (S[ispin][1] - spinProjectedOnQ[1]) * phase;
										unpolarised[2] += sign * (S[ispin][2] - spinProjectedOnQ[2]) * phase;

										nuclear += sign * phase;
									}
								}
							}
						}
					}

					phase = std::exp(i * double(Nconfig % 3) * dotProduct(Q3D, delta));
					unpolarised[0] *= phase;
					unpolarised[1] *= phase;
					unpolarised[2] *= phase;
					nuclear *= phase;


					//gF_magnetic += std::abs(norm(unpolarised));
					gF_magnetic += std::abs(norm(unpolarised)) * sintheta;
					gF_nuclear = std::abs(std::conj(nuclear) * nuclear);

					nuclear = 0;
					unpolarised[0] = 0;
					unpolarised[1] = 0;
					unpolarised[2] = 0;
				}
			}

			gF_2DPlot[ir] += gF_magnetic;

			gF_plot[ir][0] += gF_magnetic;
			gF_plot[ir][1] += gF_nuclear;
			std::cout << "gF_magnetic = " << gF_magnetic << "\tgF_nuclear = " << gF_nuclear;
			std::cout << "\t\tgF_plot[" << ir << "][0] = " << gF_plot[ir][0] << "\tgF_plot[" << ir << "][1] = " << gF_plot[ir][1] << std::endl;
			gF_magnetic = 0;
			gF_nuclear = 0;
		}

		else { //integrated On a Plane
			int itheta, iphi;
			for (iphi = 0; iphi < angle; iphi++) {

				r = gF_QlimDw + (gF_QlimUp - gF_QlimDw) * ir / double(R);

				costheta = 0;
				sintheta = 1;

				sinphi = sin(pi2 * iphi / double(angle));
				cosphi = cos(pi2 * iphi / double(angle));

				Q3D[0] = r * sintheta * cosphi;
				Q3D[1] = r * sintheta * sinphi;
				Q3D[2] = r * costheta;

				Qperp3D[0] = r * sintheta * cosphi;
				Qperp3D[1] = r * sintheta * sinphi;
				Qperp3D[2] = 0;


				if (r == 0) {
					QHat3D[0] = 0;
					QHat3D[1] = 0;
					QHat3D[2] = 0;
					QperpHat3D[0] = 0;
					QperpHat3D[1] = 0;
					QperpHat3D[2] = 0;
				}
				else {
					QHat3D[0] = sintheta * cosphi;
					QHat3D[1] = sintheta * sinphi;
					QHat3D[2] = costheta;

					QperpHat3D[0] = sintheta * cosphi;
					QperpHat3D[1] = sintheta * sinphi;
					QperpHat3D[2] = 0;
				}

				for (iX = 0; iX < L; iX++) {//loop over the real lattice
					for (iY = 0; iY < L; iY++) {
						for (itri = 0; itri < 4; itri++) {

							for (ispin = 0; ispin < 4; ispin++) {

								if (ispin == 3) { sign = -1; } //the non kagome spins are out.
								else if (ispin == g_lattice[0][iX][iY][itri]) { sign = -1; }
								else { sign = +1; }
								//sign = -1 if the spin is out, else it is +1. (because the direction of the spin in S is in,
								//and then must be flipped in function of the direction of the actual spin)

								//phase du to the spin position
								phase = std::exp(i * dotProduct(Q3D, g_SpinPosition[iX][iY][itri][ispin]));

								ScalarProduct = dotProduct(QHat3D, S[ispin]);
								spinProjectedOnQ[0] = ScalarProduct * QHat3D[0];
								spinProjectedOnQ[1] = ScalarProduct * QHat3D[1];
								spinProjectedOnQ[2] = ScalarProduct * QHat3D[2];

								//the whole spin is taken into account (unpolarised neutrons)
								unpolarised[0] += sign * (S[ispin][0] - spinProjectedOnQ[0]) * phase;
								unpolarised[1] += sign * (S[ispin][1] - spinProjectedOnQ[1]) * phase;
								unpolarised[2] += sign * (S[ispin][2] - spinProjectedOnQ[1]) * phase;

								nuclear += sign * phase;

							}
						}
					}



					//gF_magnetic += std::abs(norm(unpolarised));
					gF_magnetic += std::abs(norm(unpolarised)) * sintheta;
					gF_nuclear = std::abs(std::conj(nuclear) * nuclear);

					nuclear = 0;
					unpolarised[0] = 0;
					unpolarised[1] = 0;
					unpolarised[2] = 0;

				}
			}

			gF_2DPlot[ir] += gF_magnetic;

			gF_plot[ir][0] += gF_magnetic;
			gF_plot[ir][1] += gF_nuclear;
			//std::cout << "gF_magnetic = " << gF_magnetic << "\tgF_nuclear = " << gF_nuclear;
			//std::cout << "\t\tgF_plot[" << ir << "][0] = " << gF_plot[ir][0] << "\tgF_plot[" << ir << "][1] = " << gF_plot[ir][1] << std::endl;
			gF_magnetic = 0;
			gF_nuclear = 0;
		}

	}
	std::cout << "gF_plot[" << ir << "][0] = " << gF_plot[ir][0] << "\tgF_plot[" << ir << "][1] = " << gF_plot[ir][1] << std::endl;

	gF_layer++;
	std::cout << "It took " << time(NULL) - initial_time << " seconds" << std::endl;
}

/// <summary>
/// Save the plot in FlavienPlot.txt
/// </summary>
void endFlavienPlot() {
	int ir;
	int i, imax;
	double Q[R];

	std::fstream myfile;
	std::remove("FlavienPlot.txt");
	myfile.open("FlavienPlot.txt", std::fstream::out | std::fstream::app);
	myfile << "Q (A)" << "\t" << "magnetic2D" << "\t" << "magnetic3D" << std::endl;

	double magnetic2D[R], Plot3D[R];

	for (ir = 0; ir < R; ir++) {
		Q[ir] = gF_QlimDw + (gF_QlimUp - gF_QlimDw) * ir / double(R);
		Q[ir] *= NT / NT_HoTiO;//size of experimental cell
		Q[ir] *= 0.0000000001;//in angstrom

		gF_2DPlot[ir] /= Niteration;
		gF_2DPlot[ir] /= g_Nspins; //Nspins
		gF_2DPlot[ir] /= angle;

		gF_plot[ir][0] /= Niteration;
		gF_plot[ir][0] /= g_Nspins; //Nspins
		gF_plot[ir][0] /= angle;
		gF_plot[ir][0] /= angle;

		gF_plot[ir][1] /= Niteration;
		gF_plot[ir][1] /= g_Nspins; //Nspins
		gF_plot[ir][1] /= angle;
		gF_plot[ir][1] /= angle;

		Plot3D[ir] = 0;

	}

	for (imax = 0; imax < R; imax++) {
		for (i = 0; i <= imax; i++) {
			Plot3D[imax] += gF_2DPlot[i] * Q[i];
		}
		Plot3D[imax] /= Q[imax];
	}

	for (ir = 0; ir < R; ir++) {
		//myfile << Q[ir] << "\t " << gF_plot[ir][0] << "\t" << gF_plot[ir][0] << std::endl;
		myfile << Q[ir] << "\t " << gF_plot[ir][0] << "\t" << Plot3D[ir] << std::endl;
	}
}



/// <summary>
/// Create the spheres in reciprocal space on which FlavienPlotUpdate will integrated on.
/// </summary>
void plot2DInitialisation() {
	//gF_QlimDw = 1.0 * pi2 / NT;;
	//gF_QlimUp = 2.0 * pi2 / NT;

	gF_QlimDw = 0 * pi2 / NT;;
	gF_QlimUp = 3 * pi2 / NT;

	gF_layer = 0;

	int i0, i1;
	int ir, iphi, itheta;
	double r;
	double phi, theta;

	//intialise 2D parameters
	for (ir = 0; ir < R; ir++) {
		for (iphi = 0; iphi < angle; iphi++) {
			phi = pi2 * iphi / angle;
			r = gF_QlimDw + (gF_QlimUp - gF_QlimDw) * ir / R;

			gF_Q2D[ir][iphi][0] = r * cos(phi);
			gF_Q2D[ir][iphi][1] = r * sin(phi);
			gF_Q2D[ir][iphi][2] = 0;

			if (r < 0.00001) {
				gF_QHat2D[ir][iphi][0] = 0;
				gF_QHat2D[ir][iphi][1] = 0;
				gF_QHat2D[ir][iphi][2] = 0;
			}
			else {
				gF_QHat2D[ir][iphi][0] = cos(phi);
				gF_QHat2D[ir][iphi][1] = sin(phi);
				gF_QHat2D[ir][iphi][2] = 0;
			}

		}
	}

	g_delta[0] = ReaL[2][0];
	g_delta[1] = ReaL[2][1];
	g_delta[2] = ReaL[2][2];
}

/// <summary>
/// Integrate the intensity of the neutron scattering on the sphere of |Q| = cste.
/// </summary>
void plot2DUpdate() {
	time_t initial_time = time(NULL);
	std::cout << "Flavien's plot update n" << gF_layer << " is being calculated...";

	int ir, iphi, j;
	gF_nuclear = 0;
	gF_magnetic = 0;

	int iX, iY, itri, ispin;

	double r;
	double phi, normsqrt;
	double rng;

	double Q3D[3] = { 0, 0, 0 };
	double Qperp3D[3] = { 0, 0, 0 };

	double QHat3D[3] = { 0, 0, 0 };
	double QperpHat3D[3] = { 0, 0, 0 };

	double sign;
	double ScalarProduct;
	double spinProjectedOnQ[3] = { 0, 0, 0 };
	double delta[3];
	const std::complex<double> i(0, 1);
	std::complex<double> phase = 0;
	std::complex<double> unpolarised[3] = { 0, 0, 0 };

	std::complex<double> nuclear = 0;

	for (ir = 0; ir < R; ir++) {

		for (iphi = 0; iphi < angle; iphi++) {
			

			r = gF_QlimDw + (gF_QlimUp - gF_QlimDw) * ir / double(R);
			phi = pi2 / double(angle);

			Q3D[0] = r * cos(phi) * RecL[0][0] + r * sin(phi) * RecL[1][0];
			Q3D[1] = r * cos(phi) * RecL[0][1] + r * sin(phi) * RecL[1][1];
			Q3D[2] = r * cos(phi) * RecL[0][2] + r * sin(phi) * RecL[1][2];//*/

			/*Q3D[0] = +1.0 * r * cos(phi) - sqrt(3) * r * sin(phi);
			Q3D[1] = +1.0 * r * cos(phi) + sqrt(3) * r * sin(phi);
			Q3D[2] = -2.0 * r * cos(phi);//*/
			
			normsqrt = sqrt(Q3D[0] * Q3D[0] + Q3D[1] * Q3D[1] + Q3D[2] * Q3D[2]);

			if (normsqrt != 0) {
				QHat3D[0] = Q3D[0] / normsqrt;
				QHat3D[1] = Q3D[1] / normsqrt;
				QHat3D[2] = Q3D[2] / normsqrt;
			}
			else {
				QHat3D[0] = 0;
				QHat3D[1] = 0;
				QHat3D[2] = 0;

			}
			//std::cout << "test 1" << std::endl;

			for (iX = 0; iX < L; iX++) {//loop over the real lattice
				for (iY = 0; iY < L; iY++) {
					for (itri = 0; itri < 4; itri++) {
						for (ispin = 0; ispin < 3; ispin++) {

							if (ispin == 3) { sign = -1; } //the non kagome spins are out.
							else if (ispin == g_lattice[0][iX][iY][itri]) { sign = -1; }
							else { sign = +1; }
							//sign = -1 if the spin is out, else it is +1. (because the direction of the spin in S is in,
							//and then must be flipped as a function of the direction of the actual spin)

							//phase du to the spin position
							phase = std::exp(i * dotProduct(gF_Q2D[ir][iphi], g_SpinPosition[iX][iY][itri][ispin]));

							ScalarProduct = dotProduct(gF_QHat2D[ir][iphi], S[ispin]);
							spinProjectedOnQ[0] = ScalarProduct * QHat3D[0];
							spinProjectedOnQ[1] = ScalarProduct * QHat3D[1];
							spinProjectedOnQ[2] = ScalarProduct * QHat3D[2];

							//the whole spin is taken into account (unpolarised neutrons)
							unpolarised[0] += sign * (S[ispin][0] - spinProjectedOnQ[0]) * phase;
							unpolarised[1] += sign * (S[ispin][1] - spinProjectedOnQ[1]) * phase;
							unpolarised[2] += sign * (S[ispin][2] - spinProjectedOnQ[2]) * phase;

							nuclear += sign * phase;

						}
					}
				}
			}

			//gF_magnetic += std::abs(norm(unpolarised));-+++++++++++
			gF_2DPlot[ir] += std::abs(norm(unpolarised));

			nuclear = 0;
			unpolarised[0] = 0;
			unpolarised[1] = 0;
			unpolarised[2] = 0;
		}
	}

	gF_layer++;
	std::cout << "It took " << time(NULL) - initial_time << " seconds" << std::endl;
}

/// <summary>
/// Save the plot in FlavienPlot.txt
/// </summary>
void end2DPlot() {
	int ir;
	int i, imax;
	double Q[R];

	std::fstream myfile;
	std::remove("FlavienPlot.txt");
	myfile.open("FlavienPlot.txt", std::fstream::out | std::fstream::app);
	myfile << "Q (A)" << "\t" << "magnetic2D" << "\t" << "magnetic3D" << std::endl;

	double magnetic2D[R], Plot3D[R];

	for (ir = 0; ir < R; ir++) {
		Q[ir] = gF_QlimDw + (gF_QlimUp - gF_QlimDw) * ir / double(R);
		Q[ir] *= NT / NT_HoTiO;//size of experimental cell
		Q[ir] *= 0.0000000001;//in angstrom

		gF_2DPlot[ir] /= Niteration;
		gF_2DPlot[ir] /= g_Nspins; //Nspins
		gF_2DPlot[ir] /= angle;

		gF_plot[ir][0] /= Niteration;
		gF_plot[ir][0] /= g_Nspins; //Nspins
		gF_plot[ir][0] /= angle;

		Plot3D[ir] = 0;

	}

	for (imax = 0; imax < R; imax++) {
		for (i = 0; i <= imax; i++) {
			Plot3D[imax] += gF_2DPlot[i] * Q[i];
		}
		Plot3D[imax] /= Q[imax];
	}

	for (ir = 0; ir < R; ir++) {
		//myfile << Q[ir] << "\t " << gF_plot[ir][0] << "\t" << gF_plot[ir][0] << std::endl;
		myfile << Q[ir] << "\t" << gF_2DPlot[ir] << std::endl;
	}
}



/// <summary>
/// Calculate the order parameter between the star phase and the columnar phase.
/// </summary>
double calculateOrderParameter() {
	time_t time_stamp_1 = time(NULL);

	int ik, iX, iY, itri, ispin;
	double sign;

	const std::complex<double> i(0, 1);
	std::complex<double> phase = 0;
	std::complex<double> orderP[3] = { 0,0,0 };

	double OrderParameter = 0;

	//direction in the brilloin zone
	double Q[3] = { 2 * pi2 / (3.0 * NT), 0, 0 };

	double Qhat[3] = { 0,0,0 };
	Qhat[0] = Q[0] / (sqrt(norm(Q)));
	Qhat[1] = Q[1] / (sqrt(norm(Q)));
	Qhat[2] = Q[2] / (sqrt(norm(Q)));

	double ScalarProduct;
	std::complex<double> spinProjectedOnQ[3];

	for (ik = 0; ik < M; ik++) {
		for (iX = 0; iX < L; iX++) {//loop over the real lattice
			for (iY = 0; iY < L; iY++) {
				for (itri = 0; itri < 4; itri++) {
					for (ispin = 0; ispin < 3; ispin++) {
						if (g_lattice[ik][iX][iY][itri] == ispin) { sign = -1; }
						else { sign = +1; }

						phase = std::exp(i * dotProduct(Q, g_SpinPosition[iX][iY][itri][ispin]));
						
						ScalarProduct = dotProduct(Qhat, S[ispin]);
						spinProjectedOnQ[0] = ScalarProduct * Qhat[0];
						spinProjectedOnQ[1] = ScalarProduct * Qhat[1];

						//orderP[0] += sign * (S[ispin][0] - spinProjectedOnQ[0]) * phase;
						orderP[1] += sign * (S[ispin][1] - spinProjectedOnQ[1]) * phase;
					}
				}
			}
		}//end of the loop on the real lattice
	}

	OrderParameter = sqrt(std::abs(norm(orderP)));
	OrderParameter /= double(g_Nspins); //Turn OrderParameter from an extensive to an intensive variable.
	OrderParameter /= Sperp; //Turn OrderParameter from an extensive to an intensive variable.
	OrderParameter = 9.0 * OrderParameter / 4.0; //Turn OrderParameter is between 0 and 1.
	return(OrderParameter);
}



/// <summary>
/// Return the dot product of two vectors of length 3 made of doubles and const doubles.
/// </summary>
/// <param name="A"></param>
/// <param name="B"></param>
/// <returns></returns>
double dotProduct(double A[3], const double B[3]) {
	return(A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
}

/// <summary>
/// Return the dot product of two vectors of length 3 made of doubles.
/// </summary>
/// <param name="A"></param>
/// <param name="B"></param>
/// <returns></returns>
double dotProduct(double A[3], double B[3]) {
	return(A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
}

/// <summary>
/// Return the dot product of two vectors of length 3 made of complexes.
/// </summary>
/// <param name="A"></param>
/// <param name="B"></param>
/// <returns></returns>
std::complex<double> dotProduct(std::complex<double> A[3], std::complex<double> B[3]) {
	return(std::conj(A[0]) * B[0] + std::conj(A[1]) * B[1] + std::conj(A[2]) * B[2]);
}

/// <summary>
/// Return A' . A
/// </summary>
/// <param name="A"></param>
/// <param name="B"></param>
/// <returns></returns>
double norm(double A[3]) {
	return(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
}

/// <summary>
/// Return A' . A
/// </summary>
/// <param name="A"></param>
/// <param name="B"></param>
/// <returns></returns>
std::complex<double> norm(std::complex<double> A[3]) {
	return(std::real(std::conj(A[0]) * A[0] + std::conj(A[1]) * A[1] + std::conj(A[2]) * A[2]));
}
