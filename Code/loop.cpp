#include <iostream>
#include <stdlib.h>
#include <cmath>

#include "node.h"
#include "loop.h"
#include "defs.h"

extern int compteur;
extern int N_external_loop;
extern int g_Nspins;

loop::loop(int k) {
	m_k = k;
	m_tail = NULL;
	m_head = NULL;
}

loop::loop(int k, node* tail) {
	m_k = k;
	m_tail = tail;
	m_head = m_tail;
}

loop::loop(int k, int x, int y, char tri, int dimer)
{
	m_k = k;
	m_tail = new node(x, y, tri, dimer);
	m_head = m_tail;
}

node* loop::getHead() {
	return(m_head);
}

void loop::addNode(node* newNode) {
	if (m_tail == NULL) { m_head = newNode; m_tail = m_head; }
	else { m_head = (*m_head).addNext(newNode); }//m_head->next=newNode; m_head=newNode;
}

void loop::updateLoop(int x, int y, char tri, int dimer) {
	node* newNode = new node(x, y, tri, dimer);
	if (m_tail == NULL) { m_tail = newNode; m_head = newNode; }

	else {
		if (!(*this).updateNode(newNode)) {
			(*this).addNode(newNode);
		}
	}
}

node* loop::cutTheTail()
{
	node* newTail = (*m_tail).nextNode(); //newTail points at the second metamere.
	delete m_tail; //m_tail points at NULL and newTail points at now the first metamere.
	m_tail = newTail;
	return(m_tail);
}

bool loop::updateNode(node* newNode) {
	node* anyNode = m_tail;
	while (anyNode != NULL) {

		if ((*anyNode).isEqualTo(newNode)) {
			(*anyNode).updateDimer(newNode);
			delete newNode;
			return(true);
		}
		anyNode = (*anyNode).nextNode();
	}

	return(false);
}

void loop::loopTheWorm()
{
	while (!(*m_head).isEqualTo(m_tail)) {
		(*this).cutTheTail();
	}
	(*m_head).updateDimer(m_tail);
	(*this).cutTheTail();
}

void loop::deleteUnchangedNode(int Lattice[L][L][4]) {
	node* anyNode = m_tail;
	if (anyNode == NULL) { std::cout << "(*this) is NULL" << std::endl; }
	while (anyNode != m_head) {
		if (anyNode == m_tail) {
			if ((*anyNode).hasChangedLattice(Lattice)) {
				if ((*(*anyNode).nextNode()).hasChangedLattice(Lattice)) {
					anyNode = (*anyNode).nextNode();
				}
				else {
					(*anyNode).deleteNext();
				}
			}
			else {//delete the node (it has been corrected)
				m_tail = (*this).cutTheTail();
			}
		}
		else {
			if ((*(*anyNode).nextNode()).hasChangedLattice(Lattice)) {
				anyNode = (*anyNode).nextNode();
			}
			else {//delete the node (it has been corrected)
				if ((*anyNode).nextNode() == m_head) {
					(*anyNode).deleteNext();
					m_head = anyNode;
				}
				else { (*anyNode).deleteNext(); }
			}
		}
	}
}

void loop::showLoop()
{
	std::cout << std::endl << "showLoop" << std::endl;
	node* anyNode = m_tail;
	if (anyNode == NULL) { std::cout << "Loop is empty" << std::endl; }
	else {
		std::cout << "head : ";
		m_head->showNode();
		std::cout << std::endl << "tail : ";
		m_tail->showNode();
		while (anyNode != NULL) {
			std::cout << std::endl;
			(*anyNode).showNode();
			anyNode = (*anyNode).nextNode();
		}
	}
	std::cout << std::endl << std::endl;
}

bool loop::isNotEmpty()
{
	if (m_head == NULL || m_tail == NULL) { return(false); }
	else { return(true); }
}

bool loop::contains(int xNext, int yNext, char triNext, int dimerNext) {
	node* anyNode = m_tail;
	while (anyNode != NULL) {
		if ((*anyNode).isEqualTo(xNext, yNext, triNext, dimerNext)) { return(true); }
		anyNode = (*anyNode).nextNode();
	}
	return(false);
}

node* loop::moveNode(node* Node) {
	node* tempNode = Node->nextNode();
	this->addNode(tempNode);
	Node->updateNext(tempNode->nextNode());
	tempNode->updateNext(NULL);
	m_head = tempNode;
	return(m_head);
}

loop* loop::returnLoop(int tempLattice[L][L][4])
{
	// The initial intuition might to take the first node (m_tail) of (*this) loop then to find a neighbor in the
	// following, realocate it, and then to do it again and again until (*this) loop is empty.
	// There is a major difficulty here because of loops that are close from each other .
	// In deed, there will be some case where the two neighbors of a dimer will flipped in a single loop.
	// That why there is the function "testCorrectNeighbor" that verify if or not following the selected neighbor
	// interfere with an other flipped loop.

	loop* realLoop = new loop(this->m_k, m_tail->copy());

	node* anyNode; //anyNode is in (*this) loop

	node* neighbor1 = new node(0, 0, 'A', 0);
	node* neighbor2 = new node(0, 0, 'A', 0);

	(*m_tail).getNeighbor(neighbor1, neighbor2);

	anyNode = m_tail;


	while ((*anyNode).nextNode() != NULL) {

		if ((*(*anyNode).nextNode()).isEqualTo(neighbor1)) {

			if ((*this).testCorrectNeighbor((*anyNode).nextNode(), (*realLoop).getHead(), tempLattice)) {
				(*(*realLoop).moveNode(anyNode)).getNeighbor(neighbor1, neighbor2);
				anyNode = m_tail;
			}
			else {
				anyNode = (*anyNode).nextNode();
			}
		}

		else if ((*(*anyNode).nextNode()).isEqualTo(neighbor2)) {

			if ((*this).testCorrectNeighbor((*anyNode).nextNode(), (*realLoop).getHead(), tempLattice)) {
				(*(*realLoop).moveNode(anyNode)).getNeighbor(neighbor1, neighbor2);
				anyNode = m_tail;
			}
			else {
				anyNode = (*anyNode).nextNode();
			}
		}

		else {
			anyNode = (*anyNode).nextNode();
		}
	}


	(*this).cutTheTail();

	delete neighbor1;
	delete neighbor2;


	return(realLoop);
}

bool loop::testCorrectNeighbor(node* neighbor, node* lastAdded, int tempLattice[L][L][4])
{
	/*Both neighbor1 and neighbor2 are in (*this) loop. Which one is the next after node ?
	d stand for dimer

					   / \				   /d\
					  /   \				  /   \   (  if b1 didn't flipped, then Node's neighbor is neighbor2)
					 /Node \			 /     \  (else b2 didn't flipped,  and Node's neighbor is neighbor1)
					(In Loop)			neighbor1
				   /________d\_________/b1_______\_________
							  \  down /			  \       /
							   \ tri /			   \     /
								\   /				\   /
								 \ /				 \ /
								  b2
								 / \
								/   \
							   /     \
							  neighbor2
							 /d________\
	*/

	int b = (*lastAdded).spinInBetween(neighbor);//b and a must be different. If there're equal, there is defect in down tri

	if ((*neighbor).isBFlipped(b, m_k, tempLattice)) { return(true); }

	return(false);
}

loop::~loop()
{
	while (m_tail != NULL) { (*this).cutTheTail(); }
}