#include <iostream>

#include "defs.h"
#include "node.h"
#include "next.h"


node::node(int x, int y, char tri, int dimer) : m_x(x), m_y(y), m_tri(tri), m_dimer(dimer)
{
	m_next = NULL;
}

node* node::nextNode()
{
	return(m_next);
}

node* node::addNext(node* newNode)
{
	m_next = newNode;
	return(newNode);
}

void node::deleteNext() {
	node* tempNode = m_next;
	m_next = tempNode->m_next;
	//std::cout << "delete is called.";
	delete tempNode;
}

void node::updateNext(node* Node) {
	m_next = Node;
}

bool node::isEqualTo(node* M)
{
	if ((m_x == M->m_x) && (m_y == M->m_y) && (m_tri == M->m_tri)) { return(true); } //doesn't concern the dimer !
	else { return(false); }
}

bool node::isEqualTo(int xNext, int yNext, char triNext, int dimerNext)
{
	if ((m_x == xNext) && (m_y == yNext) && (m_tri == triNext) && (m_dimer == dimerNext)) { return(true); } //doesn't concern the dimer !
	return(false);
}

bool node::sameDimer(node* M)
{
	if (M->m_dimer == m_dimer) { return(true); }
	else { return(false); }
}

void node::updateDimer(node* newDimer)
{
	m_dimer = newDimer->m_dimer;
}

void node::getNeighbor(node* neighbor1, node* neighbor2)
{
	if (m_dimer == 0) {
		next(m_x, m_y, m_tri, m_dimer, 1, 0, &(neighbor1->m_x), &(neighbor1->m_y), &(neighbor1->m_tri), &(neighbor1->m_dimer));
		next(m_x, m_y, m_tri, m_dimer, 2, 0, &(neighbor2->m_x), &(neighbor2->m_y), &(neighbor2->m_tri), &(neighbor2->m_dimer));
	}
	else if (m_dimer == 1) {
		next(m_x, m_y, m_tri, m_dimer, 0, 0, &(neighbor1->m_x), &(neighbor1->m_y), &(neighbor1->m_tri), &(neighbor1->m_dimer));
		next(m_x, m_y, m_tri, m_dimer, 2, 0, &(neighbor2->m_x), &(neighbor2->m_y), &(neighbor2->m_tri), &(neighbor2->m_dimer));
	}
	else {
		next(m_x, m_y, m_tri, m_dimer, 0, 0, &(neighbor1->m_x), &(neighbor1->m_y), &(neighbor1->m_tri), &(neighbor1->m_dimer));
		next(m_x, m_y, m_tri, m_dimer, 1, 0, &(neighbor2->m_x), &(neighbor2->m_y), &(neighbor2->m_tri), &(neighbor2->m_dimer));
	}
}

node* node::copy()
{
	return(new node(m_x, m_y, m_tri, m_dimer));
}

bool node::hasChangedLattice(int Lattice[L][L][4])
{
	switch (m_tri)
	{
	case 'A':
		if (Lattice[m_x][m_y][0] == m_dimer) { return(false); }
		return(true);
	case 'B':
		if (Lattice[m_x][m_y][1] == m_dimer) { return(false); }
		return(true);
	case 'C':
		if (Lattice[m_x][m_y][2] == m_dimer) { return(false); }
		return(true);
	case 'D':
		if (Lattice[m_x][m_y][3] == m_dimer) { return(false); }
		return(true);
	}
	return(false);
}

void node::showNode() {
	std::cout << m_x << " " << m_y << " " << m_tri << " " << m_dimer << " " << m_next;
}

int node::spinInBetween(node* neighbor) {
	return(findSpinBetweenTwoOthers(m_tri, m_dimer, neighbor->m_tri));
}

int node::getDimer() {
	return(m_dimer);
}

void node::findNextSpin(int c, int nnew, int* xout, int* yout, char* triout, int* dimerout) {
	next(m_x, m_y, m_tri, c, nnew, 0, xout, yout, triout, dimerout);
}

bool node::isBFlipped(int b, int k, int tempLattice[L][L][4])
{
	switch (m_tri)
	{
	case 'A':
		return(isItFlipped(m_x, m_y, 0, b, k, tempLattice));
	case 'B':
		return(isItFlipped(m_x, m_y, 1, b, k, tempLattice));
	case 'C':
		return(isItFlipped(m_x, m_y, 2, b, k, tempLattice));
	case 'D':
		return(isItFlipped(m_x, m_y, 3, b, k, tempLattice));
	}
	return(false);
}

node::~node() {
	//std::cout << "node deleted" << std::endl;
}
