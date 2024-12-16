#include <iostream>

#include "check.h"

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

bool check(int x, int y, char tri, int n) {
	if (x == xA && y == yA && tri == triA && n == nA) {
		return(true);
	}
	if (x == xB && y == yB && tri == triB && n == nB) {
		return(true);
	}
	if (x == xC && y == yC && tri == triC && n == nC) {
		return(true);
	}
	return(false);
}