#ifndef DEFS_H
#define DEFS_H

//-------------ALGORITHM PARAMETERS-------------
//Defines the lattice size in the real space
constexpr int L=18;

//Defines the lattice size in the Trotter dimension
constexpr int M = 1; //500;

//Defines the magnetic field// Field norm
constexpr double Bi = 0.3533;
constexpr double Phi = 0; // 30° = 0.52359877559; // in rad

//Defines the interaction constant
constexpr double g = 0;
constexpr double Vi = -1;
constexpr double delta = -0;

//Defines the temperature.
constexpr double Ti = 0.005; //Initial temperature

//Defines the number of MC steps needed to equilibrate the lattice.
constexpr int Neq = 2000; // Equilibration done before each latticeMeasure() call.




//-------------MEASURES PARAMETERS-------------
//Defines the temperature domain.
constexpr double Tf = 1; //Final temperature measured
constexpr double dT = 1; //Temperature step

//Defines the temperature domain.
constexpr double Bf = 0.3538; //Final magnetic field measured
constexpr double dB = 0.00001; //Magnetic fiel step

//Defines the temperature domain.
constexpr double Vf = -2; //Final diagonal term measured
constexpr double dV = -10; //Diagonal term step

//Defines the number of configuration use to measure the system equilibrium properties.
constexpr int Nmeas = 10000;

//There is approximatly Nmeas*Neq loops proposed during the algorithm

//--------------Zacharie code------------------
constexpr bool DoYouWantZacharieCode = false;
//Defines if we take the random reference to initialise the lattice
constexpr bool DoYouWantLatticeRef = true;
//Defines if we save the current final lattice as a reference
constexpr bool DoYouWantToSaveRef = false;

constexpr bool DoYouWantTransitionGraphNSP = false;
constexpr int Nconfiguration = 5;


//--------------Magnetisation probability distribution PARAMETERS-------------------
constexpr bool magnetizationHistogram = false;

//Defines the number of random draw of magnetisation.
constexpr int Nrd = 10000;



//--------------NSP PARAMETERS-------------------
constexpr bool DoYouWantNeutronScatteringPlot = false; //true = yes    false = no

//Defines the number of NSP measures.
constexpr int NSpinConfig = 20;

//defines NSP Q's mesh.
constexpr int xmult = 2;
constexpr int ymult = 2;
constexpr int xdiv = xmult * 1 * L; 
constexpr int ydiv = ymult * 1 * L;


//------------2D Plot---------------------
constexpr bool DoYouWant2DPlot = false;

//------------Flavien's Plot---------------------
constexpr bool DoYouWantFlavienPlot = false;
constexpr int Ndirection = 100;
constexpr int angle = 96;
constexpr int R = 250; //Number of points in the Flavien's plot.
constexpr int Niteration = 50;
constexpr double NT_HoTiO = 0.000000000714;
constexpr bool ferromagneticPeak = false;
constexpr bool integratedOnASphere = true;



//DO NOT TOUCH :( (Lattice constants)
constexpr double S[4][3] = { { 0 , -0.94280904158 , 0.33333333333 }, { +0.81649658092 , 0.47140452079 , 0.33333333333 }, { -0.81649658092 , 0.47140452079 , 0.33333333333  }, {0, 0, -1} };//Kagome spins
constexpr double pi2 = 2 * 3.14159265358979323846; //2*pi

constexpr double Sperp = 0.94280904158; //Sperp = 2 * sqrt(2) / 3

constexpr double NT = 1 / (1.41421356237); //distance between two up triangles, NT means Nearest Triange NT = 1 / sqrt(2);

constexpr double ReaL[3][3] = { {    NT    ,            0                , 0 }, { NT / 2 ,      NT * 1.73205080757 / 2    , 0  }, { NT / 2, NT / (2 * 1.73205080757), NT * 0.81649658092} };//ReaL lattice basis     0.81649658092 = sqrt(2/3)
constexpr double RecL[3][3] = { { pi2 / NT , -pi2 / (NT * 1.73205080757) , -pi2 / (NT * 2.44948974278) }, {   0    , 2 * pi2 / (NT * 1.73205080757) , -pi2 / (NT * 2.44948974278) }, { 0, 0, pi2 * 1.73205080757 / (NT * 1.41421356237)} };//Reciprocal lattice basis     2.44948974278 = sqrt(6)

#endif //DEFS_H
