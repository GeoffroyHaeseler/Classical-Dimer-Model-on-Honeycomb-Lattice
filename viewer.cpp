#include <iostream>
#include <fstream> 
#include <math.h>
#include <string>
#include <sstream> 

#include "defs.h"
#include "viewer.h"

extern int g_lattice[L][L][L][4];

void view() {
	viewLattice();
	int k;
	for (k = 0; k < M; k++) { viewResult(k); }
}



void viewLattice() {
	std::fstream myfile;
	remove("lattice.ply");
	myfile.open("lattice.ply", std::fstream::out | std::fstream::app);

	myfile << "ply" << std::endl;
	myfile << "format ascii 1.0" << std::endl;
	myfile << "element vertex " << 24 * L * L << std::endl;
	myfile << "property float x" << std::endl;
	myfile << "property float y" << std::endl;
	myfile << "property float z" << std::endl;
	myfile << "element face " << 8 * L * L << std::endl;
	myfile << "property list uchar int vertex_index" << std::endl;
	myfile << "end_header" << std::endl;
	myfile.close();

	int ix, iy;
	for (ix = 0; ix < L; ix++) {
		for (iy = 0; iy < L; iy++) {
			addLatticeVertices(ix, iy);
		}
	}

	for (ix = 0; ix < L; ix++) {
		for (iy = 0; iy < L; iy++) {
			addLatticeFaces(ix, iy);
		}
	}

}

void addLatticeVertices(int ix, int iy) {
	double r3 = sqrt(3);
	double r3_2 = 0.5 * sqrt(3);
	double X, Y;

	std::fstream myfile;
	myfile.open("lattice.ply", std::fstream::out | std::fstream::app);

	//Triangle A up and down
	X = 4.0 * ix;
	Y = 4.0 * r3_2 * iy;

	myfile << X + 0.5 << " " << Y + r3_2 << " " << 0 << "    {Vertex " << int(ix) << " " << int(iy) << " A}" << std::endl;
	myfile << X << " " << Y << " " << 0 << std::endl;
	myfile << X + 1 << " " << Y << " " << 0 << std::endl;

	myfile << X + 0.5 << " " << Y + r3_2 << " " << 0 << std::endl;
	myfile << X + 1 << " " << Y + r3 << " " << 0 << std::endl;
	myfile << X << " " << Y + r3 << " " << 0 << std::endl;

	//Triangle B up and down
	X = 4.0 * ix + 2;
	Y = 4.0 * r3_2 * iy;

	myfile << X + 0.5 << " " << Y + r3_2 << " " << 0 << "    {Vertex " << int(ix) << " " << int(iy) << " B}" << std::endl;
	myfile << X << " " << Y << " " << 0 << std::endl;
	myfile << X + 1 << " " << Y << " " << 0 << std::endl;

	myfile << X + 0.5 << " " << Y + r3_2 << " " << 0 << std::endl;
	myfile << X + 1 << " " << Y + r3 << " " << 0 << std::endl;
	myfile << X << " " << Y + r3 << " " << 0 << std::endl;

	//Triangle C up and down
	X = 4.0 * ix + 3.0;
	Y = 4.0 * r3_2 * iy + r3;

	myfile << X + 0.5 << " " << Y + r3_2 << " " << 0 << "    {Vertex " << int(ix) << " " << int(iy) << " C}" << std::endl;
	myfile << X << " " << Y << " " << 0 << std::endl;
	myfile << X + 1 << " " << Y << " " << 0 << std::endl;

	myfile << X + 0.5 << " " << Y + r3_2 << " " << 0 << std::endl;
	myfile << X + 1 << " " << Y + r3 << " " << 0 << std::endl;
	myfile << X << " " << Y + r3 << " " << 0 << std::endl;

	//Triangle D up and down
	X = 4.0 * ix + 1.0;
	Y = 4.0 * r3_2 * iy + r3;

	myfile << X + 0.5 << " " << Y + r3_2 << " " << 0 << "    {Vertex " << int(ix) << " " << int(iy) << " D}" << std::endl;
	myfile << X << " " << Y << " " << 0 << std::endl;
	myfile << X + 1 << " " << Y << " " << 0 << std::endl;

	myfile << X + 0.5 << " " << Y + r3_2 << " " << 0 << std::endl;
	myfile << X + 1 << " " << Y + r3 << " " << 0 << std::endl;
	myfile << X << " " << Y + r3 << " " << 0 << std::endl;

	myfile.close();
}

void addLatticeFaces(int ix, int iy) {
	std::fstream myfile;
	myfile.open("lattice.ply", std::fstream::out | std::fstream::app);

	int N = 24 * (iy + ix * L);

	//Only now N is is number of the line {ix iy iz tetrahedron + S0}

	myfile << 3 << " " << N + 0 << " " << N + 1 << " " << N + 2 << "    {Face " << ix << " " << iy << " A}" << std::endl;
	myfile << 3 << " " << N + 3 << " " << N + 4 << " " << N + 5 << std::endl;
	myfile << 3 << " " << N + 6 << " " << N + 7 << " " << N + 8 << std::endl;
	myfile << 3 << " " << N + 9 << " " << N + 10 << " " << N + 11 << std::endl;
	myfile << 3 << " " << N + 12 << " " << N + 13 << " " << N + 14 << std::endl;
	myfile << 3 << " " << N + 15 << " " << N + 16 << " " << N + 17 << std::endl;
	myfile << 3 << " " << N + 18 << " " << N + 19 << " " << N + 20 << std::endl;
	myfile << 3 << " " << N + 21 << " " << N + 22 << " " << N + 23 << std::endl;

	myfile.close();
}



void viewResult(int k) {

	int ix, iy;

	std::stringstream ss;
	ss << "results_" << k << ".ply";
	std::string str = ss.str();
	const char* name = str.c_str();

	std::fstream myfile;
	remove(name);
	myfile.open(name, std::fstream::out | std::fstream::app);

	myfile << "ply" << std::endl;
	myfile << "format ascii 1.0" << std::endl;
	myfile << "element vertex " << 36 * L * L << std::endl;
	myfile << "property float x" << std::endl;
	myfile << "property float y" << std::endl;
	myfile << "property float z" << std::endl;
	myfile << "element face " << 12 * L * L << std::endl;
	myfile << "property list uchar int vertex_index" << std::endl;
	myfile << "end_header" << std::endl;
	myfile.close();

	for (ix = 0; ix < L; ix++) {
		for (iy = 0; iy < L; iy++) {
			//std::cout << std::endl << k << " A" << g_lattice[k][ix][iy][0] << " B" << g_lattice[k][ix][iy][1] << " C" << g_lattice[k][ix][iy][2] << " D" << g_lattice[k][ix][iy][3] << std::endl;
			addArrowsVertices(k, ix, iy, name);
		}
	}

	for (ix = 0; ix < L; ix++) {
		for (iy = 0; iy < L; iy++) {
			addArrowsFaces(k, ix, iy, name);
		}
	}
}

void addArrowsVertices(int k, int ix, int iy, const char* name) {
	std::fstream myfile;
	myfile.open(name, std::fstream::out | std::fstream::app);

	double r3 = sqrt(3);
	double r3_2 = 0.5 * sqrt(3);
	double r3_10 = 0.1 * sqrt(3);
	double r3_20 = 0.05 * sqrt(3);
	double X, Y;
	float sign;

	//Triangle A spin 0
	X = 4 * ix + 0.5;
	Y = 4 * r3_2 * iy + r3_2;
	sign = 1;
	if (g_lattice[k][ix][iy][0] == 1 || g_lattice[k][ix][iy][0] == 2) { sign = -1; }
	myfile << X << " " << Y + sign * 0.2 << " " << 0 << "    {Spin " << int(ix) << " " << int(iy) << " A 0}" << std::endl;
	myfile << X - sign * 0.1 << " " << Y - sign * 0.2 << " " << 0 << std::endl;
	myfile << X + sign * 0.1 << " " << Y - sign * 0.2 << " " << 0 << std::endl;

	//Triangle A spin 1
	X = 4.0 * ix;
	Y = 4.0 * r3_2 * iy;
	sign = 1;
	if (g_lattice[k][ix][iy][0] == 1) { sign = -1; }
	myfile << X + sign * r3_10 << " " << Y + sign * 0.1 << " " << 0 << "    {Spin " << int(ix) << " " << int(iy) << " A 1}" << std::endl;
	myfile << X - sign * (0.05 + r3_10) << " " << Y + sign * (r3_20 - 0.1) << " " << 0 << std::endl;
	myfile << X + sign * (0.05 - r3_10) << " " << Y - sign * (r3_20 + 0.1) << " " << 0 << std::endl;

	//Triangle A spin 2
	X = 4.0 * ix + 1.0;
	Y = 4.0 * r3_2 * iy;
	sign = 1;
	if (g_lattice[k][ix][iy][0] == 2) { sign = -1; }
	myfile << X - sign * r3_10 << " " << Y + sign * 0.1 << " " << 0 << "    {Spin " << int(ix) << " " << int(iy) << " A 1}" << std::endl;
	myfile << X - sign * (0.05 - r3_10) << " " << Y - sign * (r3_20 + 0.1) << " " << 0 << std::endl;
	myfile << X + sign * (0.05 + r3_10) << " " << Y + sign * (r3_20 - 0.1) << " " << 0 << std::endl;



	//Triangle B spin 0
	X = 4.0 * ix + 2.0 + 0.5;
	Y = 4.0 * r3_2 * iy + r3_2;
	sign = 1;
	if (g_lattice[k][ix][iy][1] == 1 || g_lattice[k][ix][iy][1] == 2) { sign = -1; }
	myfile << X << " " << Y + sign * 0.2 << " " << 0 << "    {Spin " << int(ix) << " " << int(iy) << " A 0}" << std::endl;
	myfile << X - sign * 0.1 << " " << Y - sign * 0.2 << " " << 0 << std::endl;
	myfile << X + sign * 0.1 << " " << Y - sign * 0.2 << " " << 0 << std::endl;

	//Triangle B spin 1
	X = 4.0 * ix + 2.0;
	Y = 4.0 * r3_2 * iy;
	sign = 1;
	if (g_lattice[k][ix][iy][1] == 1) { sign = -1; }
	myfile << X + sign * r3_10 << " " << Y + sign * 0.1 << " " << 0 << "    {Spin " << int(ix) << " " << int(iy) << " A 1}" << std::endl;
	myfile << X - sign * (0.05 + r3_10) << " " << Y + sign * (r3_20 - 0.1) << " " << 0 << std::endl;
	myfile << X + sign * (0.05 - r3_10) << " " << Y - sign * (r3_20 + 0.1) << " " << 0 << std::endl;

	//Triangle B spin 2
	X = 4.0 * ix + 2.0 + 1.0;
	Y = 4.0 * r3_2 * iy;
	sign = 1;
	if (g_lattice[k][ix][iy][1] == 2) { sign = -1; }
	myfile << X - sign * r3_10 << " " << Y + sign * 0.1 << " " << 0 << "    {Spin " << int(ix) << " " << int(iy) << " A 1}" << std::endl;
	myfile << X - sign * (0.05 - r3_10) << " " << Y - sign * (r3_20 + 0.1) << " " << 0 << std::endl;
	myfile << X + sign * (0.05 + r3_10) << " " << Y + sign * (r3_20 - 0.1) << " " << 0 << std::endl;



	//Triangle C spin 0
	X = 4.0 * ix + 3.0 + 0.5;
	Y = 4.0 * r3_2 * iy + r3 + r3_2;
	sign = 1;
	if (g_lattice[k][ix][iy][2] == 1 || g_lattice[k][ix][iy][2] == 2) { sign = -1; }
	myfile << X << " " << Y + sign * 0.2 << " " << 0 << "    {Spin " << int(ix) << " " << int(iy) << " A 0}" << std::endl;
	myfile << X - sign * 0.1 << " " << Y - sign * 0.2 << " " << 0 << std::endl;
	myfile << X + sign * 0.1 << " " << Y - sign * 0.2 << " " << 0 << std::endl;

	//Triangle C spin 1
	X = 4.0 * ix + 3.0;
	Y = 4.0 * r3_2 * iy + r3;
	sign = 1;
	if (g_lattice[k][ix][iy][2] == 1) { sign = -1; }
	myfile << X + sign * r3_10 << " " << Y + sign * 0.1 << " " << 0 << "    {Spin " << int(ix) << " " << int(iy) << " A 1}" << std::endl;
	myfile << X - sign * (0.05 + r3_10) << " " << Y + sign * (r3_20 - 0.1) << " " << 0 << std::endl;
	myfile << X + sign * (0.05 - r3_10) << " " << Y - sign * (r3_20 + 0.1) << " " << 0 << std::endl;

	//Triangle C spin 2
	X = 4.0 * ix + 3.0 + 1.0;
	Y = 4.0 * r3_2 * iy + r3;
	sign = 1;
	if (g_lattice[k][ix][iy][2] == 2) { sign = -1; }
	myfile << X - sign * r3_10 << " " << Y + sign * 0.1 << " " << 0 << "    {Spin " << int(ix) << " " << int(iy) << " A 1}" << std::endl;
	myfile << X - sign * (0.05 - r3_10) << " " << Y - sign * (r3_20 + 0.1) << " " << 0 << std::endl;
	myfile << X + sign * (0.05 + r3_10) << " " << Y + sign * (r3_20 - 0.1) << " " << 0 << std::endl;



	//Triangle D spin 0
	X = 4.0 * ix + 1.0 + 0.5;
	Y = 4.0 * r3_2 * iy + r3 + r3_2;
	sign = 1;
	if (g_lattice[k][ix][iy][3] == 1 || g_lattice[k][ix][iy][3] == 2) { sign = -1; }
	myfile << X << " " << Y + sign * 0.2 << " " << 0 << "    {Spin " << int(ix) << " " << int(iy) << " A 0}" << std::endl;
	myfile << X - sign * 0.1 << " " << Y - sign * 0.2 << " " << 0 << std::endl;
	myfile << X + sign * 0.1 << " " << Y - sign * 0.2 << " " << 0 << std::endl;

	//Triangle D spin 1
	X = 4.0 * ix + 1.0;
	Y = 4.0 * r3_2 * iy + r3;
	sign = 1;
	if (g_lattice[k][ix][iy][3] == 1) { sign = -1; }
	myfile << X + sign * r3_10 << " " << Y + sign * 0.1 << " " << 0 << "    {Spin " << int(ix) << " " << int(iy) << " A 1}" << std::endl;
	myfile << X - sign * (0.05 + r3_10) << " " << Y + sign * (r3_20 - 0.1) << " " << 0 << std::endl;
	myfile << X + sign * (0.05 - r3_10) << " " << Y - sign * (r3_20 + 0.1) << " " << 0 << std::endl;

	//Triangle D spin 2
	X = 4.0 * ix + 1.0 + 1.0;
	Y = 4.0 * r3_2 * iy + r3;
	sign = 1;
	if (g_lattice[k][ix][iy][3] == 2) { sign = -1; }
	myfile << X - sign * r3_10 << " " << Y + sign * 0.1 << " " << 0 << "    {Spin " << int(ix) << " " << int(iy) << " A 1}" << std::endl;
	myfile << X - sign * (0.05 - r3_10) << " " << Y - sign * (r3_20 + 0.1) << " " << 0 << std::endl;
	myfile << X + sign * (0.05 + r3_10) << " " << Y + sign * (r3_20 - 0.1) << " " << 0 << std::endl;

	myfile.close();
}

void addArrowsFaces(int k, int ix, int iy, const char* name) {
	std::fstream myfile;
	myfile.open(name, std::fstream::out | std::fstream::app);

	int N = 36 * (iy + ix * L);

	//Only now N is is number of the line {ix iy iz tetrahedron + S0}

	myfile << 3 << " " << N + 0 << " " << N + 1 << " " << N + 2 << "    {Arrow " << ix << " " << iy << " A 0}" << std::endl;
	myfile << 3 << " " << N + 3 << " " << N + 4 << " " << N + 5 << std::endl;
	myfile << 3 << " " << N + 6 << " " << N + 7 << " " << N + 8 << std::endl;
	myfile << 3 << " " << N + 9 << " " << N + 10 << " " << N + 11 << std::endl;
	myfile << 3 << " " << N + 12 << " " << N + 13 << " " << N + 14 << std::endl;
	myfile << 3 << " " << N + 15 << " " << N + 16 << " " << N + 17 << std::endl;
	myfile << 3 << " " << N + 18 << " " << N + 19 << " " << N + 20 << std::endl;
	myfile << 3 << " " << N + 21 << " " << N + 22 << " " << N + 23 << std::endl;
	myfile << 3 << " " << N + 24 << " " << N + 25 << " " << N + 26 << std::endl;
	myfile << 3 << " " << N + 27 << " " << N + 28 << " " << N + 29 << std::endl;
	myfile << 3 << " " << N + 30 << " " << N + 31 << " " << N + 32 << std::endl;
	myfile << 3 << " " << N + 33 << " " << N + 34 << " " << N + 35 << std::endl;

	myfile.close();
}

