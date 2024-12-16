void next(const int xin, const int yin, const char triin, const int nin, const int nnew, const int local, int* xout, int* yout, char* triout, int* nout);

void next(const int xin, const int yin, const int triin, const int nin, const int nnew, const int local, int* xout, int* yout, int* triout, int* nout);

/// <summary>
/// Return the spin between the dimer 1 and the dimer 2. The dimer 1 points in the shared down triangle.
/// </summary>
/// <param name="x1"></param>
/// <param name="y1"></param>
/// <param name="tri1"></param>
/// <param name="n1"></param>
/// <param name="x2"></param>
/// <param name="y2"></param>
/// <param name="tri2"></param>
/// <param name="n2"></param>
/// <returns></returns>
int findSpinBetweenTwoOthers(const char tri1, const int n1, const char tri2);

/// <summary>
/// Return true if the dimer in x, y and tri is flipped between g_lattice[k] and tempLattice. If not return false.
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="tri"></param>
/// <param name="c"></param>
/// <param name="k"></param>
/// <param name="tempLattice"></param>
/// <returns></returns>
bool isItFlipped(int x, int y, int tri, int dimer, int k, int tempLattice[L][L][4]);
