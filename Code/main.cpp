#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <fstream>


#include "defs.h"
#include "initialise.h"
#include "addAWorm.h"
#include "qfFunctions.h"
#include "measure.h"
#include "next.h"
#include "viewer.h"

extern int g_lattice[M][L][L][4];
extern int g_latticeTransitionGraphLattice[L][L][4][3];

int main() {
	time_t initial_time = time(NULL);
	time_t time_stamp_1;
	initialise();

	double T, B, V;
	int i, j;

	/*int ik, ix, iy, ispin;
	for (ix = 0; ix < L; ix++) {
	for (iy = 0; iy < L; iy++) {
		for (ispin = 0; ispin < 4; ispin++) {
			std::cout << g_lattice[0][ix][iy][ispin] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}*/



	//Neutron Scattering Plot
	if (DoYouWantNeutronScatteringPlot) {
		probaInitialisation(Ti, Bi, Vi, Phi);//set probas for selected parameters
		for (i = 0; i < NSpinConfig; i++) {
			for (j = 0; j < 50 * Neq; j++) { addAQuantumWorm(); }
			latticeMeasure(Vi, 0);
		}
		endMeasure(Ti, Bi, Vi, Phi, Nmeas);
	}

	else if (DoYouWantZacharieCode) {

		if (DoYouWantTransitionGraphNSP) {
			probaInitialisation(Ti, Bi, Vi, Phi);//set probas for selected parameters
			int ix, iy, itri;

			for (i = 0; i < Nconfiguration; i++) {

				for (j = 0; j < Neq; j++) { addAQuantumWorm(); }

				for (ix = 0; ix < L; ix++) {
					for (iy = 0; iy < L; iy++) {
						for (itri = 0; itri < 4; itri++) {

							g_latticeTransitionGraphLattice[ix][iy][itri][(g_lattice[0][ix][iy][itri] + 0) % 3] = -1;
							g_latticeTransitionGraphLattice[ix][iy][itri][(g_lattice[0][ix][iy][itri] + 1) % 3] = +1;
							g_latticeTransitionGraphLattice[ix][iy][itri][(g_lattice[0][ix][iy][itri] + 2) % 3] = +1;

						}
					}
				}

				for (j = 0; j < Neq; j++) { addAQuantumWorm(); }

				for (ix = 0; ix < L; ix++) {
					for (iy = 0; iy < L; iy++) {
						for (itri = 0; itri < 4; itri++) {

							g_latticeTransitionGraphLattice[ix][iy][itri][(g_lattice[0][ix][iy][itri] + 0) % 3] -= -1;
							g_latticeTransitionGraphLattice[ix][iy][itri][(g_lattice[0][ix][iy][itri] + 1) % 3] -= +1;
							g_latticeTransitionGraphLattice[ix][iy][itri][(g_lattice[0][ix][iy][itri] + 2) % 3] -= +1;

						}
					}
				}

				latticeMeasure(Vi, 0);
			}
			endMeasure(Ti, Bi, Vi, Phi, Nmeas);
		}

		else {
			probaInitialisation(Ti, Bi, Vi, Phi);//set probas for selected parameters

			for (j = 0; j < 1000; j++) { addAQuantumWorm(); }

			std::cout << std::endl;
			std::cout << std::endl;

			save_lattice();

		}
	}

	//Neutron Scattering of 3d powder
	else if (DoYouWantFlavienPlot) {
		probaInitialisation(Ti, Bi, Vi, Phi);//set probas for selected parameters
		for (i = 0; i < Niteration; i++) {
			for (j = 0; j < 100 * Neq; j++) { addAQuantumWorm(); }
			latticeMeasure(Vi, i);
		}
		endMeasure(Ti, Bi, Vi, Phi, Nmeas);
	}

	//Neutron Scattering of 2d powder
	else if (DoYouWant2DPlot) {
		probaInitialisation(Ti, Bi, Vi, Phi);//set probas for selected parameters

		for (i = 0; i < Niteration; i++) {
			time_stamp_1 = time(NULL);
			for (j = 0; j < 100 * Neq; j++) { addAQuantumWorm(); }
			latticeMeasure(Vi, i);
			std::cout << "The loop " << i << " took " << time(NULL) - time_stamp_1 << " seconds to complete" << std::endl << std::endl << std::endl;
		}
		endMeasure(Ti, Bi, Vi, Phi, Nmeas);
	}

	//Magnetisation probability distribution
	else if (magnetizationHistogram) {
		probaInitialisation(Ti, Bi, Vi, Phi);//set probas for selected parameters
		for (j = 0; j < 100000; j++) { addAQuantumWorm(); }
		std::cout << "test";
		while(true) {
			for (j = 0; j < Neq; j++) { addAQuantumWorm(); }
			orderParameterHistogramMeasure();
		}
	}

	//Thermodynamical study
	else {
	
		for (T = Ti; T <= Tf; T += dT) {
			for (B = Bi; B <= Bf; B += dB) {
				for (V = Vi; V >= Vf; V += dV) {
					time_stamp_1 = time(NULL);
					probaInitialisation(T, B, V, Phi);//set probas for selected parameters

					for (j = 0; j < 1000 * Neq; j++) { addAQuantumWorm(); }

					//Lattice measure
					for (i = 0; i < Nmeas; i++) {

						for (j = 0; j < Neq; j++) { addAQuantumWorm(); }

						latticeMeasure(V, 0);
					} //do the measures

					endMeasure(T, B, V, Phi, Nmeas);

					std::cout << "Measurments took " << time(NULL) - time_stamp_1 << " seconds to complete" << std::endl;
					std::cout << std::endl;
					std::cout << std::endl;
				}
			}
		}



		/*for (T = 0.75; T < 1.5; T += 0.025) {
			for (B = Bi; B <= Bf; B += dB) {
				for (V = Vi; V >= Vf; V += dV) {
					time_stamp_1 = time(NULL);
					probaInitialisation(T, B, V, Phi);//set probas for selected parameters

					for (j = 0; j < 10 * Neq; j++) { addAQuantumWorm(); }

					//Lattice measure
					for (i = 0; i < Nmeas; i++) {

						for (j = 0; j < Neq; j++) { addAQuantumWorm(); }

						latticeMeasure(V, 0);
					} //do the measures

					endMeasure(T, B, V, Phi, Nmeas);

					std::cout << "Measurments took " << time(NULL) - time_stamp_1 << " seconds to complete" << std::endl;
					std::cout << std::endl;
					std::cout << std::endl;
				}
			}
		}

		for (T = 1.5; T < 2; T += 0.05) {
			for (B = Bi; B <= Bf; B += dB) {
				for (V = Vi; V >= Vf; V += dV) {
					time_stamp_1 = time(NULL);
					probaInitialisation(T, B, V, Phi);//set probas for selected parameters

					for (j = 0; j < 10 * Neq; j++) { addAQuantumWorm(); }

					//Lattice measure
					for (i = 0; i < Nmeas; i++) {

						for (j = 0; j < Neq; j++) { addAQuantumWorm(); }

						latticeMeasure(V, 0);
					} //do the measures

					endMeasure(T, B, V, Phi, Nmeas);

					std::cout << "Measurments took " << time(NULL) - time_stamp_1 << " seconds to complete" << std::endl;
					std::cout << std::endl;
					std::cout << std::endl;
				}
			}
		}

		for (T = 2; T <= 5; T += 0.2) {
			for (B = Bi; B <= Bf; B += dB) {
				for (V = Vi; V >= Vf; V += dV) {
					time_stamp_1 = time(NULL);
					probaInitialisation(T, B, V, Phi);//set probas for selected parameters

					for (j = 0; j < 10 * Neq; j++) { addAQuantumWorm(); }

					//Lattice measure
					for (i = 0; i < Nmeas; i++) {

						for (j = 0; j < Neq; j++) { addAQuantumWorm(); }

						latticeMeasure(V, 0);
					} //do the measures

					endMeasure(T, B, V, Phi, Nmeas);

					std::cout << "Measurments took " << time(NULL) - time_stamp_1 << " seconds to complete" << std::endl;
					std::cout << std::endl;
					std::cout << std::endl;
				}
			}
		}*/
	
	}


	std::cout << "The program took " << time(NULL) - initial_time << " seconds to complete" << std::endl << std::endl << std::endl;
	return(0);
}
