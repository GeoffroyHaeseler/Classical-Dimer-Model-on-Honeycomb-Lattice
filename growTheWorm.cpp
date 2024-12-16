#include <iostream>     /*in/out console print*/
#include <stdlib.h>     /* srand, rand */
#include <cmath>

#include "growTheWorm.h"
#include "defs.h"
#include "next.h"

extern double z0;
extern double z1;
extern double z2;
extern double kT;
extern double proba[3][3][3][3][3][3][3];
extern int g_lattice[M][L][L][4];
extern bool g_TypeAHexagone[L][L][4];
extern bool g_TypeCHexagone[L][L][4];

/// <summary>
/// Choose the next triangle. If it is a worm piece, return false : we need to close the loop. If not, return true
/// </summary>
/// <param name="Worm1s"></param>
/// <returns></returns>
void growTheWorm(int k, int* x, int* y, char* tri, int* n) {
	/*


		  Configuration from addAWorm: the bottom defect (o) is in a down triangle


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
		 /_________\         /_________\




		 First we move the defect to the up triangle from which the dimer come from
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
		   \     /             \     /
			\   /               \   /
			 \ /                 \ /
			 /d\                 / \
			/   \               /   \
		   /     \             /  o  \
		  /       \           /       \
		 /_________\         /_________\


		 The we randomly move the defect, here right:
				   A_________B          _________
					\       /           \       /
					 \  o  /             \  o  /
					  \   /               \   /
					   \ /                 \ /
					   / \                 /d\
					  /   \               /   \
					 /     \             /     \
					/       \           /       \
		  _________/________d\_________/_________\
		  \       /           \       /
		   \     /             \     /
			\   /               \   /
			 \ /                 \ /
			 /d\                 / \
			/   \               /   \
		   /     \             /     \
		  /       \           /       \
		 /_________\         /________d\


	*/

	//std::cout << "down tri : x=" << *x << "\ty=" << *y << "\ttri=" << *tri << "\tn=" << *n << std::endl;
	//x, y, z, tri and n indicates the position and the new direction of the last dimer flip
	int prev_n = *n; //contains the previous direction

	//We need to pass the down triangle to go into the next up triangle. To do this, we must follow the dimer.
	passTheDownTriangle(k, x, y, tri, n);


	//std::cout << "up tri : x=" << *x << "\ty=" << *y << "\ttri=" << *tri << "\tn=" << *n << "       ";

	//We must now choose which down triangle Worm1 will enter.
	double p0, p1, p2;
	calculateProbability(k, *x, *y, *tri, &p0, &p1, &p2);


	double rng = (double)rand() / RAND_MAX;
	//We choose spin 0
	if (0 <= rng && rng < p0) { next(*x, *y, *tri, *n, 0, 1, x, y, tri, n); }

	//We choose spin 1
	else if (p0 <= rng && rng < p0 + p1) { next(*x, *y, *tri, *n, 1, 1, x, y, tri, n); }

	//We choose spin 2
	else { next(*x, *y, *tri, *n, 2, 1, x, y, tri, n); }

	//std::cout << "We left the triangle through " << *n << std::endl;

	//we are now able to identify a new node. We still didn't add it to the worm which stay in the previous down triangle.
	switch (*tri)
	{
	case 'A':
		g_lattice[k][*x][*y][0] = *n;
		break;
	case 'B':
		g_lattice[k][*x][*y][1] = *n;
		break;
	case 'C':
		g_lattice[k][*x][*y][2] = *n;
		break;
	case 'D':
		g_lattice[k][*x][*y][3] = *n;
		break;
	}
}

/// <summary>
/// (k, x, y, tri, n) point at the last flipped dimer. At the end of this function, it will point at the next dimer that must be flipped.
/// </summary>
/// <param name="k"></param>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="tri"></param>
/// <param name="n"></param>
void passTheDownTriangle(int k, int* x, int* y, char* tri, int* n) {
	bool isDownTriPassed = false;
	int look_x, look_y, look_dimer, i;
	char look_tri;
	//We are in a down triangle, we reach the next up triangle while following the 2 out 1 in rule
	for (i = 0; i < 3; i++) {//We look the three spins arround us. Each one is coming from a different up triangle.
		if (*n != i && !isDownTriPassed) {//If it's not the one we come from
			next(*x, *y, *tri, *n, i, 0, &look_x, &look_y, &look_tri, &look_dimer); //we take up triangle's coordinates,
			//	     next dimer __/    \__ local = 0 : we change triangle.
			switch (look_tri)
			{
			case 'A':
				if (g_lattice[k][look_x][look_y][0] == look_dimer) {//and if it is the down triangle dimer
					*x = look_x;//we leave the down triangle through it.
					*y = look_y;
					*tri = look_tri;
					*n = i;
					//std::cout << look_x << " " << look_y << " " << look_tri << i << std::endl;
					isDownTriPassed = true;
				}
				break;
			case 'B':
				if (g_lattice[k][look_x][look_y][1] == look_dimer) {//and if it is the down triangle dimer
					*x = look_x;//we leave the down triangle through it.
					*y = look_y;
					*tri = look_tri;
					*n = i;
					//std::cout << look_x << " " << look_y << " " << look_tri << i << std::endl;
					isDownTriPassed = true;
				}
				break;
			case 'C':
				if (g_lattice[k][look_x][look_y][2] == look_dimer) {//and if it is the down triangle dimer
					*x = look_x;//we leave the down triangle through it.
					*y = look_y;
					*tri = look_tri;
					*n = i;
					//std::cout << look_x << " " << look_y << " " << look_tri << i << std::endl;
					isDownTriPassed = true;
				}
				break;
			case 'D':
				if (g_lattice[k][look_x][look_y][3] == look_dimer) {//and if it is the down triangle dimer
					*x = look_x;//we leave the down triangle through it.
					*y = look_y;
					*tri = look_tri;
					*n = i;
					//std::cout << look_x << " " << look_y << " " << look_tri << i << std::endl;
					isDownTriPassed = true;
				}
				break;
			}
		}
	}

}

/// <summary>
/// Calculate the probability for the monopole to go out of the (k, x, y, tri) triangle, knowing it enters through g_lattice[k][x][y][tri]
/// </summary>
/// <param name="k: label the magnet we're in in the Trotter dimension."></param>
/// <param name="x: label the cell position we're in."></param>
/// <param name="y: label the cell position we're in."></param>
/// <param name="tri: label the triangle we're in."></param>
/// <param name="p0: output the propaibility to go to 0."></param>
/// <param name="p1: output the propaibility to go to 1."></param>
/// <param name="p2: output the propaibility to go to 2."></param>
void calculateProbability(const int k, int x, int y, char tri, double* p0, double* p1, double* p2)
{
	/*
	0 is given by x, y and tri = triNum
	   /B\___/C\
	  /   \ /   \
	/A\___/0\___/D\
	   \ /   \ /
	   /F\___/E\
	*/

	int triNum = 0;
	switch (tri)
	{
	case'A':
		triNum = 0;
		break;
	case'B':
		triNum = 1;
		break;
	case'C':
		triNum = 2;
		break;
	case'D':
		triNum = 3;
		break;
	}

	//CALCULATE THE IMPACT OF THE PLAQUETTE TERM

	int Vcount0    = 0;
	int deltacount0 = 0;
	int Vcount1    = 0;
	int deltacount1 = 0;
	int Vcount2    = 0;
	int deltacount2 = 0;

	int N; //doesn't mean anything
	//neighbour positions
	int XA, YA, TriA;
	int XB, YB, TriB;
	int XC, YC, TriC;
	int XD, YD, TriD;
	int XE, YE, TriE;
	int XF, YF, TriF;

	//positions of the six triangles around
	next(x, y, triNum, 1, 2, 0, &XA, &YA, &TriA, &N);//A
	next(x, y, triNum, 0, 2, 0, &XB, &YB, &TriB, &N);//B
	next(x, y, triNum, 0, 1, 0, &XC, &YC, &TriC, &N);//C
	next(x, y, triNum, 2, 1, 0, &XD, &YD, &TriD, &N);//D
	next(x, y, triNum, 2, 0, 0, &XE, &YE, &TriE, &N);//E
	next(x, y, triNum, 1, 0, 0, &XF, &YF, &TriF, &N);//F

	//calculate the diagonal energy (V) du to the dimer being in
	//0
	if (g_lattice[k][XA][YA][TriA] == 2 && g_lattice[k][XB][YB][TriB] == 1) {
		Vcount0++;
		if (g_TypeAHexagone[XA][YA][TriA]) { deltacount0--; }
		if (g_TypeCHexagone[XA][YA][TriA]) { deltacount0++; }
	}
	if (g_lattice[k][XC][YC][TriC] == 2 && g_lattice[k][XD][YD][TriD] == 1) {
		Vcount0++;
		if (g_TypeAHexagone[x ][y ][tri ]) { deltacount0--; }
		if (g_TypeCHexagone[x ][y ][tri ]) { deltacount0++; }
	}

	//1
	if (g_lattice[k][XE][YE][TriE] == 0 && g_lattice[k][XF][YF][TriF] == 2) {
		Vcount1++;
		if (g_TypeAHexagone[XF][YF][TriF]) { deltacount1--; }
		if (g_TypeCHexagone[XF][YF][TriE]) { deltacount1++; }
	}
	if (g_lattice[k][XA][YA][TriA] == 0 && g_lattice[k][XB][YB][TriB] == 2) {
		Vcount1++;
		if (g_TypeAHexagone[XA][YA][TriA]) { deltacount1--; }
		if (g_TypeCHexagone[XA][YA][TriA]) { deltacount1++; }
	}

	//2
	if (g_lattice[k][XC][YC][TriC] == 1 && g_lattice[k][XD][YD][TriD] == 0) {
		Vcount2++;
		if (g_TypeAHexagone[x ][y ][tri ]) { deltacount2--; }
		if (g_TypeCHexagone[x ][y ][tri ]) { deltacount2++; }
	}
	if (g_lattice[k][XE][YE][TriE] == 1 && g_lattice[k][XF][YF][TriF] == 0) {
		Vcount2++;
		if (g_TypeAHexagone[XF][YF][TriF]) { deltacount2--; }
		if (g_TypeCHexagone[XF][YF][TriE]) { deltacount2++; }
	}	
	deltacount0++; deltacount1++; deltacount2++; //turn them into a integer between 0 and 2 rather than between -1 and 1 to use them as array 


	//update the probabilities
	*p0 = proba[Vcount0][Vcount1][Vcount2][deltacount0][deltacount1][deltacount2][0];	//P(0->0)
	*p1 = proba[Vcount0][Vcount1][Vcount2][deltacount0][deltacount1][deltacount2][1];	//P(0->1)
	*p2 = proba[Vcount0][Vcount1][Vcount2][deltacount0][deltacount1][deltacount2][2];	//P(0->2)



	//std::cout << "p0=" << *p0 << "\tp1=" << *p1 << "\tp2=" << *p2 << std::endl;
	//those probabilities can be above 1 or below 0 if V!=0, it's not an issue, since their sum is equal to one.
}
