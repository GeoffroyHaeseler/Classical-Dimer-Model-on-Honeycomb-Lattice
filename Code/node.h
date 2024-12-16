#ifndef DEF_NODE
#define DEF_NODE

#include "defs.h"
// The metamere is a single segment of the worm body.
// with this analogy, the worm is made of a list of metameres.
// (it has neither head nor tail : they just are respctivly the last and the first element of the list)
// Here, each metamere is an up tetrahedron recognised by it's position (m_x, m_y, m_z, m_tetra),
// it dimer direction (m_dimer) and the next metamere in the worm (m_next).

class node
{

public:

	/// <summary>
	/// Create a node with the specified parameters
	/// </summary>
	/// <param name="x"></param> x coordinate of the cell
	/// <param name="y"></param> y coordinate of the cell
	/// <param name="tetra"></param> tetrahedron of the metamere
	/// <param name="dimer"></param> direction of the dimer
	node(int x, int y, char tri, int dimer); //overloading

	/// <summary>
	/// returning the next node.
	/// </summary>
	/// <returns></returns>
	node* nextNode();

	/// <summary>
	/// Return the value of the next node
	/// </summary>
	/// <param name="newMetamere"></param>
	node* addNext(node* newNode);

	/// <summary>
	/// Delete the next Node if not NULL. Also link the Nodes correctly.
	/// </summary>
	/// <param name="tempNode"></param>
	void deleteNext();

	/// <summary>
	/// this->m_next = Node. Do not deallocate memory.
	/// </summary>
	void updateNext(node* Node);

	/// <summary>
	/// Send back true if the parameter and the indifier indicates the same triangle. If not, send back false.
	/// </summary>
	/// <param name="M"></param>
	/// <returns></returns>
	bool isEqualTo(node* M);

	/// <summary>
	/// Return true if (*this) node contains *Next dimer at the exact same position). Else return false.
	/// </summary>
	/// <param name="xNext"></param>
	/// <param name="yNext"></param>
	/// <param name="triNext"></param>
	/// <param name="dimerNext"></param>
	/// <returns></returns>
	bool isEqualTo(int xNext, int yNext, char triNext, int dimerNext);

	/// <summary>
	/// return true if (*this) and M have the same dimer. else return false
	/// </summary>
	/// <param name="M"></param>
	/// <returns></returns>
	bool sameDimer(node* M);

	/// <summary>
	/// Update a dimer that the loop already goes into
	/// </summary>
	/// <param name="M"></param>
	void updateDimer(node* newDimer);

	/// <summary>
	/// Create the same node but allocated elsewhere and without m_next
	/// </summary>
	/// <param name="x"></param> x coordinate of the cell
	/// <param name="y"></param> y coordinate of the cell
	/// <param name="tetra"></param> tetrahedron of the metamere
	/// <param name="dimer"></param> direction of the dimer
	node* copy();

	/// <summary>
	/// Does this node change the lattice?
	/// </summary>
	bool hasChangedLattice(int Lattice[L][L][4]);

	/// <summary>
	/// Show the node in the console.
	/// </summary>
	void showNode();

	/// <summary>
	/// Create two nodes corresponding to neighbors of (this*)
	/// </summary>
	void getNeighbor(node* neighbor1, node* neighbor2);

	/// <summary>
	/// Return the spin between (*this) and Neighbor nodes. The dimer of (*this) point toward the down triangle between the nodes.
	/// </summary>
	/// <param name="Neighbor"></param>
	/// <returns></returns>
	int spinInBetween(node* neighbor);

	/// <summary>
	/// Return the dimer of (*this) node
	/// </summary>
	/// <returns></returns>
	int getDimer();

	/// <summary>
	/// Put location of the spin next to (*this) node in the pointers.
	/// </summary>
	void findNextSpin(int c, int nnew, int* xout, int* yout, char* triout, int* dimerout);

	/// <summary>
	/// Return true if the dimer c in this -> m_x, m_y, m_tri is flipped (between g_lattice[k] and tempLattice). If not return false.
	/// </summary>
	/// <param name="c"></param>
	/// <param name="m_k"></param>
	/// <param name="tempLattice"></param>
	/// <returns></returns>
	bool isBFlipped(int c, int k, int tempLattice[L][L][4]);


	~node();

private:
	int m_x;
	int m_y;
	char m_tri;
	int m_dimer;
	node* m_next;
};
#endif
