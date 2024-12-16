#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>

#include "addAWorm.h"
#include "defs.h"
#include "next.h"
#include "check.h"
#include "loop.h"
#include "growTheWorm.h"


extern int g_lattice[M][L][L][4];
extern double g_proba[3][3];

extern int wormSize;

extern int xA;
extern int yA;
extern char triA;
extern int nA;

extern int xB;
extern int yB;
extern char triB;
extern int nB;

extern int xC;
extern int yC;
extern char triC;
extern int nC;

/// <summary>
/// Add a loop of flipped spins. Might measure the system.
/// </summary>
/// <param name="doWeMeasure"></param> If true, measure the thermodynamic variables, else, don't
loop* addAWorm(int k) {
	/*
	At the begining, the two defects are in the same triangle.
					_________           _________
					\       /           \       /
					 \     /             \     /
					  \   /               \   /
					   \ /                 \ /
					   /d\                 /d\
					  /   \               /   \
					 / o o \             /     \
					/       \           /       \
		  _________/_________\_________/_________\
		  \       /           \       /
		   \     /             \     /
			\   /               \   /
			 \ /                 \ /
			 /d\                 /d\
			/   \               /   \
		   /     \             /     \
		  /       \           /       \
		 /_________\         /_________\


	//A dimer change direction : the two defects are separate. One won't move anymore, the other will shuffle the lattice.
	//The immobile defect position is saved in the 'initial defect' parameters together with xa ya tria na xb yb trib and nb.
	//The parameters points toward the two spins that can end the loop if they are flip later on
				   A_________B          _________
					\       /           \       /
					 \  o  /             \     /
					  \   /               \   /
					   \ /                 \ /
					   / \                 /d\
					  /   \               /   \
					 /     \             /     \
					/       \           /       \
		  _________/________d\_________/_________\
		  \       /           \       /
		   \     /             \  o  /
			\   /               \   /
			 \ /                 \ /
			 /d\                 /d\
			/   \               /   \
		   /     \             /     \
		  /       \           /       \
		 /_________A         B_________\



	Adam's thesis notation:
			 /o\
			/   \
		   /     \
		  /       \
		 L_________R

	Present code notation:
			 /0\
			/   \
		   /     \
		  /       \
		 1_________2



	*/

	//int length = 0;
	//First, we place the head of the worm somewhere on the lattice.
	//First the coordinates of the cell
	int x = rand() % L;
	int y = rand() % L;

	//And then the triangle
	int rng = rand() % 4;
	char tri = 'A';
	switch (rng)
	{
	case 0:
		tri = 'A';
		break;
	case 1:
		tri = 'B';
		break;
	case 2:
		tri = 'C';
		break;
	case 3:
		tri = 'D';
		break;
	}
	int iniX = x;
	int iniY = y;
	char initri = tri;

	//std::cout << "We choose x=" << x << "  y=" << y <<	"  tri=" << tri << std::endl;

	//We can now identify the three position that merge the loop
	int initialDefect = g_lattice[k][x][y][rng];

	if      (initialDefect == 0) { nA = 1; nB = 2; nC = 0; }

	else if (initialDefect == 1) { nA = 0; nB = 2; nC = 1; }

	else if (initialDefect == 2) { nA = 0; nB = 1; nC = 2; }

	//We identify the three spins that are on the same down triangle as the first to use for closing check. one of these must be the last in the loop
	//std::cout << "defect is at : " << initialDefect << "   " << nA << " " << nB << " " << std::endl;
	next(x, y, tri, initialDefect, nA, 0, &xA, &yA, &triA, &nA);
	next(x, y, tri, initialDefect, nB, 0, &xB, &yB, &triB, &nB);
	xC = x;
	yC = y;
	triC = tri;

	//length++;

	//We need to choose where the spin will exit
	double p0, p1, p2;
	double rng2 = (double)rand() / RAND_MAX;
	int dimer = 0;
	calculateProbability(k, x, y, tri, &p0, &p1, &p2);
	//We choose spin 0
	if (0 <= rng2 && rng2 < p0) { dimer = 0; }

	//We choose spin 1
	else if (p0 <= rng2 && rng2 < p0 + p1) { dimer = 1; }

	//We choose spin 2
	else if (p0 + p1 <= rng2 && rng2 <= p0 + p1 + p2) { dimer = 2; }

	//We can update the lattice.
	g_lattice[k][x][y][rng] = dimer;

	//passTheDownTriangle(k, &x, &y, &tri, &dimer);


	while (!check(x, y, tri, dimer)) {
		growTheWorm(k, &x, &y, &tri, &dimer);
		//(*Loop).updateLoop(x, y, tri, dimer); //Add the node if it wasn't already changed, if it was, 1/ if it return to initial value, erase it 2/ if it is a change
		//length++;
	}
	//std::cout << length << std::endl;
	return(NULL);
}
