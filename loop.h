#ifndef DEF_LOOP
#define DEF_LOOP

#include "node.h"
#include "defs.h"
// The loop is made of a list of node.
// in our case each node is an up triangle.

class loop
{

public:
	/// <summary>
	/// Default constructor
	/// </summary>
	/// <param name="k"></param>
	loop(int k);

	/// <summary>
	/// Construct a loop in the plane k with a single node tail
	/// </summary>
	loop(int k, node* tail);

	/// <summary>
	/// Create an empty worm.
	/// </summary>
	loop(int k, int x, int y, char tri, int dimer); //constructor

	/// <summary>
	/// Return the m_head addresse
	/// </summary>
	/// <returns></returns>
	node* getHead();

	/// <summary>
	/// Create a metamere with the specified parameters and add it to the list.
	/// </summary>
	/// <param name="x"></param> x coordinate of the cell
	/// <param name="y"></param> y coordinate of the cell
	/// <param name="tetra"></param> triangle in the cell
	/// <param name="dimer"></param> direction of the dimer
	void addNode(node* newNode);

	/// <summary>
	/// Update the loop with the rule
	/// A/ Add the node if it wasn't already changed
	/// B/ if it was, update the previous dimer
	/// </summary>
	/// <param name="x"></param>
	/// <param name="y"></param>
	/// <param name="tri"></param>
	/// <param name="dimer"></param>
	void updateLoop(int x, int y, char tri, int dimer);

	/// <summary>
	/// Destroy the first metamere (corresponding to the tail).
	/// </summary>
	/// <param name="x"></param>
	/// <param name="y"></param>
	/// <param name="z"></param>
	/// <param name="tetra"></param>
	/// <param name="dimer"></param>
	node* cutTheTail();

	/// <summary>
	/// Return false if newNode is not in the loop, if newNode is in the loop modify it with the new dimer and return true.
	/// </summary>
	/// <param name="newNode"></param>
	/// <returns></returns>
	bool updateNode(node* newNode);

	/// <summary>
	/// Destroy the metameres from the tail to the first metamere that have the same properties than the head.
	/// </summary>
	/// <param name="x"></param>
	/// <param name="y"></param>
	/// <param name="z"></param>
	/// <param name="tetra"></param>
	/// <param name="dimer"></param>
	void loopTheWorm();

	/// <summary>
	/// Delete all (this) nodes that have not changed the lattice
	/// </summary>
	/// <param name="Lattice"></param>
	void deleteUnchangedNode(int Lattice[L][L][4]);

	/// <summary>
	/// Show the worm in the console
	/// </summary>
	void showLoop();

	/// <summary>
	/// Return false if (*this) is not empty. Else return true.
	/// </summary>
	/// <returns></returns>
	bool isNotEmpty();

	/// <summary>
	/// Return true if Node is in (*this) loop. Else, return false.
	/// </summary>
	/// <param name="Node"></param>
	/// <returns></returns>
	bool contains(int xout, int yout, char triout, int dimerout);

	/// <summary>
	/// this.addNode(Node->m_next); (Node->m_next) = (Node->m_next) ->m_next
	/// </summary>
	/// <param name="Node"></param>
	node* moveNode(node* Node);

	/// <summary>
	/// Return false if (*this) is not empty. Else return true.
	/// </summary>
	/// <returns></returns>
	loop* returnLoop(int tempLattice[L][L][4]);

	/// <summary>
	/// Return true if neighbor is a good candidat tobe the next Node in realLoop. Else, return false.
	/// </summary>
	/// <param name="neighbor"></param>
	/// <param name="Node"></param>
	/// <returns></returns>
	bool testCorrectNeighbor(node* neighbor, node* Node, int tempLattice[L][L][4]);

	~loop();


private:
	int m_k;
	node* m_tail;
	node* m_head;
};
#endif