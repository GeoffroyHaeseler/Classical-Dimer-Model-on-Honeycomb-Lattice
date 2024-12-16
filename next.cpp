#include <iostream>
#include <stdio.h>

#include "defs.h"
#include "next.h"

extern int g_lattice[M][L][L][4];

void next(const int xin, const int yin, const char triin, const int nin, const int nnew, const int local, int* xout, int* yout, char* triout, int* nout)
{

    /*
    xin, yin correspond to the cell coordinates in the canonical basis.
    triin correspound to a specific triangle in the cube.
    Theses information define the worm head position.
    nin correspond to a specific spin in the triangle.
    This information depend on the previous worm position.

    But we also need to know where the worm will go.
    To do this we need to know it's next location :
    The worm might stay in the same up triangle, in this case local=1, in the other case local=0.
    It will also change spin, so the user must entered the next nout, which is made through nnew.


             ___/D\___/C\___/D\___/C\
      /|\    \ /   \ /   \ /   \ /
       |     /A\___/B\___/A\___/B\___
       |        \ /   \ /   \ /   \ /
       |     ___/D\___/C\___/D\___/C\
       |     \ /   \ /   \ /   \ /
    -> |     /A\___/B\___/A\___/B\___
    ey |        \ /   \ /   \ /   \ /
       |     ___/2\___/0\___/1\___/2\
       |     \ /   \ /   \ /   \ /
       |     /0\___/1\___/2\___/0\___
       |        \ /   \ /   \ /   \ /
       |
       0---------------------------->
                    ->
                    ex

primitive cell:
             ___/D\___/C\
             \ /   \ /
             /A\___/B\___
                \ /   \ /

    */


    //The out spin will alway be the asked one.
    *nout = nnew;

    //If the triangle is not changed, nothing really change except the spin
    if (local == 1)
    {
        *xout = xin;
        *yout = yin;
        *triout = triin;
    }

    else //local == 0 => The triangle is changed
    {
        switch (triin)
        {
        case 'A':
            switch (nin) {
            case 0: //Handle the cases where we leave a triangle A from a spin 0
                switch (nnew) {
                case 1: //From triangle A to D, from a spin 0 to 1, follow -ux+uy axis
                    *xout = xin;
                    *yout = yin;
                    *triout = 'D';
                    break;
                case 2: //From triangle A to C, from a spin 0 to 2, follow  ux+uy axis
                    *xout = xin - 1;
                    *yout = yin;
                    *triout = 'C';
                    break;
                }
                break;

            case 1: //Handle the cases where we leave a triangle A from a spin 1
                switch (nnew) {
                case 0: //From triangle A to C, from a spin 1 to 0
                    *xout = xin - 1;
                    *yout = yin - 1;
                    *triout = 'C';
                    break;
                case 2: //From triangle A to B, from a spin 1 to 2
                    *xout = xin - 1;
                    *yout = yin;
                    *triout = 'B';
                    break;
                }
                break;

            case 2: //Handle the cases where we leave a triangle A from a spin 2
                switch (nnew) {
                case 0: //From triangle A to D, from a spin 2 to 0
                    *xout = xin;
                    *yout = yin - 1;
                    *triout = 'D';
                    break;
                case 1: //From triangle A to B, from a spin 2 to 1
                    *xout = xin;
                    *yout = yin;
                    *triout = 'B';
                    break;
                }
                break;
            }
            break;

        case 'B':
            switch (nin) {
            case 0: //Handle the cases where we leave a triangle B from a spin 0
                switch (nnew) {
                case 1: //From triangle B to C, from a spin 0 to 1
                    *xout = xin;
                    *yout = yin;
                    *triout = 'C';
                    break;
                case 2: //From triangle B to D, from a spin 0 to 2
                    *xout = xin;
                    *yout = yin;
                    *triout = 'D';
                    break;
                }
                break;

            case 1: //Handle the cases where we leave a triangle B from a spin 1
                switch (nnew) {
                case 0: //From triangle B to D, from a spin 1 to 0
                    *xout = xin;
                    *yout = yin - 1;
                    *triout = 'D';
                    break;
                case 2: //From triangle B to A, from a spin 1 to 2
                    *xout = xin;
                    *yout = yin;
                    *triout = 'A';
                    break;
                }
                break;

            case 2: //Handle the cases where we leave a triangle B from a spin 2
                switch (nnew) {
                case 0: //From triangle B to C, from a spin 2 to 0
                    *xout = xin;
                    *yout = yin - 1;
                    *triout = 'C';
                    break;
                case 1: //From triangle B to A, from a spin 2 to 1
                    *xout = xin + 1;
                    *yout = yin;
                    *triout = 'A';
                    break;
                }
                break;
            }
            break;

        case 'C':
            switch (nin) {
            case 0: //Handle the cases where we leave a triangle C from a spin 0
                switch (nnew) {
                case 1: //From triangle C to A, from a spin 0 to 1
                    *xout = xin + 1;
                    *yout = yin + 1;
                    *triout = 'A';
                    break;
                case 2: //From triangle C to B, from a spin 0 to 2
                    *xout = xin;
                    *yout = yin + 1;
                    *triout = 'B';
                    break;
                }
                break;

            case 1: //Handle the cases where we leave a triangle C from a spin 1
                switch (nnew) {
                case 0: //From triangle C to B, from a spin 1 to 0
                    *xout = xin;
                    *yout = yin;
                    *triout = 'B';
                    break;
                case 2: //From triangle C to D, from a spin 1 to 2
                    *xout = xin;
                    *yout = yin;
                    *triout = 'D';
                }
                break;

            case 2: //Handle the cases where we leave a triangle C from a spin 2
                switch (nnew) {
                case 0: //From triangle C to A, from a spin 2 to 0
                    *xout = xin + 1;
                    *yout = yin;
                    *triout = 'A';
                    break;
                case 1: //From triangle C to D, from a spin 2 to 1
                    *xout = xin + 1;
                    *yout = yin;
                    *triout = 'D';
                    break;
                }
                break;
            }
            break;

        case 'D':
            switch (nin) {
            case 0: //Handle the cases where we leave a triangle D from a spin 0
                switch (nnew) {
                case 1: //From triangle D to B, from a spin 0 to 1
                    *xout = xin;
                    *yout = yin + 1;
                    *triout = 'B';
                    break;
                case 2: //From triangle D to A, from a spin 0 to 2
                    *xout = xin;
                    *yout = yin + 1;
                    *triout = 'A';
                    break;
                }
                break;

            case 1: //Handle the cases where we leave a triangle D from a spin 1
                switch (nnew) {
                case 0: //From triangle D to A, from a spin 1 to 0
                    *xout = xin;
                    *yout = yin;
                    *triout = 'A';
                    break;
                case 2: //From triangle D to C, from a spin 1 to 2
                    *xout = xin - 1;
                    *yout = yin;
                    *triout = 'C';
                    break;
                }
                break;

            case 2: //Handle the cases where we leave a triangle D from a spin 2
                switch (nnew) {
                case 0: //From triangle D to B, from a spin 2 to 0
                    *xout = xin;
                    *yout = yin;
                    *triout = 'B';
                    break;
                case 1: //From triangle D to C, from a spin 2 to 1
                    *xout = xin;
                    *yout = yin;
                    *triout = 'C';
                    break;
                }
                break;
            }
            break;
        } //end switch(triin)
    }   //end else(local=0)



    //We must now take in account the PBC.
    //Since we took a complicated cell made of 4 triangles (16 spins) the PBC are trivial.
    //We just need to reset the extrem value to the other extrem :
    if (L % 2 == 0) {
        if (*yout == L) {
            *yout = 0;
            *xout -= int(L / 2);
        }
        if (*yout == -1) {
            *yout = L - 1;
            *xout += int(L / 2);
        }
    }
    else if (L % 2 == 1) {
        if (*yout == L) {
            *yout = 0;
            *xout -= int(L / 2);
            switch (*triout) {
            case 'A':
                *triout = 'B'; *xout -= 1; break;
            case 'B':
                *triout = 'A';            break;
            case 'C':
                *triout = 'D';            break;
            case 'D':
                *triout = 'C'; *xout -= 1; break;
            }
        }
        else if (*yout == -1) {
            *yout = L - 1;
            *xout += int(L / 2);
            switch (*triout) {
            case 'A':
                *triout = 'B';            break;
            case 'B':
                *triout = 'A'; *xout += 1; break;
            case 'C':
                *triout = 'D'; *xout += 1; break;
            case 'D':
                *triout = 'C';            break;
            }
        }

    }


    if (*xout >= L) { *xout -= L; }

    if (*xout < 0) { *xout += L; }


    //std::cout << xin << " " << yin << " " << triin << " becomes " << *xout << " " << *yout << " " << *triout << "\twanted: " << nin << " " << nnew << "\n";


}

void next(const int xin, const int yin, const int triin, const int nin, const int nnew, const int local, int* xout, int* yout, int* triout, int* nout)
{

    /*
    xin, yin correspond to the cell coordinates in the canonical basis.
    triin correspound to a specific triangle in the cube.
    Theses information define the worm head position.
    nin correspond to a specific spin in the triangle.
    This information depend on the previous worm position.

    But we also need to know where the worm will go.
    To do this we need it's next location :
    The worm might stay in the same up triangle, in this case local=1, in the other case local=0.
    It will also change spin, so the user must entered the next nout, which is made through nnew.


             ___/D\___/C\___/D\___/C\
      /|\    \ /   \ /   \ /   \ /
       |     /A\___/B\___/A\___/B\___
       |        \ /   \ /   \ /   \ /
       |     ___/D\___/C\___/D\___/C\
       |     \ /   \ /   \ /   \ /
    -> |     /A\___/B\___/A\___/B\___
    ey |        \ /   \ /   \ /   \ /
       |     ___/D\___/C\___/D\___/C\
       |     \ /   \ /   \ /   \ /
       |     /A\___/B\___/A\___/B\___
       |        \ /   \ /   \ /   \ /
       |
       0---------------------------->
                    ->
                    ex

primitive cell:
             ___/D\___/C\
             \ /   \ /
             /A\___/B\___
                \ /   \ /

    */


    //The out spin will alway be the asked one.
    * nout = nnew;

    //If the triangle is not changed, nothing reaaly change except the spin
    if (local == 1)
    {
        *xout = xin;
        *yout = yin;
        *triout = triin;
    }

    else //local == 0 => The triangle is changed
    {
        switch (triin)
        {
        case 0:
            switch (nin) {
            case 0: //Handle the cases where we leave a triangle A from a spin 0
                switch (nnew) {
                case 1: //From triangle A to D, from a spin 0 to 1, follow -ux+uy axis
                    *xout = xin;
                    *yout = yin;
                    *triout = 3;
                    break;
                case 2: //From triangle A to C, from a spin 0 to 2, follow  ux+uy axis
                    *xout = xin - 1;
                    *yout = yin;
                    *triout = 2;
                    break;
                }
                break;

            case 1: //Handle the cases where we leave a triangle A from a spin 1
                switch (nnew) {
                case 0: //From triangle A to C, from a spin 1 to 0
                    *xout = xin - 1;
                    *yout = yin - 1;
                    *triout = 2;
                    break;
                case 2: //From triangle A to B, from a spin 1 to 2
                    *xout = xin - 1;
                    *yout = yin;
                    *triout = 1;
                    break;
                }
                break;

            case 2: //Handle the cases where we leave a triangle A from a spin 2
                switch (nnew) {
                case 0: //From triangle A to D, from a spin 2 to 0
                    *xout = xin;
                    *yout = yin - 1;
                    *triout = 3;
                    break;
                case 1: //From triangle A to B, from a spin 2 to 1
                    *xout = xin;
                    *yout = yin;
                    *triout = 1;
                    break;
                }
                break;
            }
            break;

        case 1:
            switch (nin) {
            case 0: //Handle the cases where we leave a triangle B from a spin 0
                switch (nnew) {
                case 1: //From triangle B to C, from a spin 0 to 1
                    *xout = xin;
                    *yout = yin;
                    *triout = 2;
                    break;
                case 2: //From triangle B to D, from a spin 0 to 2
                    *xout = xin;
                    *yout = yin;
                    *triout = 3;
                    break;
                }
                break;

            case 1: //Handle the cases where we leave a triangle B from a spin 1
                switch (nnew) {
                case 0: //From triangle B to D, from a spin 1 to 0
                    *xout = xin;
                    *yout = yin - 1;
                    *triout = 3;
                    break;
                case 2: //From triangle B to A, from a spin 1 to 2
                    *xout = xin;
                    *yout = yin;
                    *triout = 0;
                    break;
                }
                break;

            case 2: //Handle the cases where we leave a triangle B from a spin 2
                switch (nnew) {
                case 0: //From triangle B to C, from a spin 2 to 0
                    *xout = xin;
                    *yout = yin - 1;
                    *triout = 2;
                    break;
                case 1: //From triangle B to A, from a spin 2 to 1
                    *xout = xin + 1;
                    *yout = yin;
                    *triout = 0;
                    break;
                }
                break;
            }
            break;

        case 2:
            switch (nin) {
            case 0: //Handle the cases where we leave a triangle C from a spin 0
                switch (nnew) {
                case 1: //From triangle C to A, from a spin 0 to 1
                    *xout = xin + 1;
                    *yout = yin + 1;
                    *triout = 0;
                    break;
                case 2: //From triangle C to B, from a spin 0 to 2
                    *xout = xin;
                    *yout = yin + 1;
                    *triout = 1;
                    break;
                }
                break;

            case 1: //Handle the cases where we leave a triangle C from a spin 1
                switch (nnew) {
                case 0: //From triangle C to B, from a spin 1 to 0
                    *xout = xin;
                    *yout = yin;
                    *triout = 1;
                    break;
                case 2: //From triangle C to D, from a spin 1 to 2
                    *xout = xin;
                    *yout = yin;
                    *triout = 3;
                }
                break;

            case 2: //Handle the cases where we leave a triangle C from a spin 2
                switch (nnew) {
                case 0: //From triangle C to A, from a spin 2 to 0
                    *xout = xin + 1;
                    *yout = yin;
                    *triout = 0;
                    break;
                case 1: //From triangle C to D, from a spin 2 to 1
                    *xout = xin + 1;
                    *yout = yin;
                    *triout = 3;
                    break;
                }
                break;
            }
            break;

        case 3:
            switch (nin) {
            case 0: //Handle the cases where we leave a triangle D from a spin 0
                switch (nnew) {
                case 1: //From triangle D to B, from a spin 0 to 1
                    *xout = xin;
                    *yout = yin + 1;
                    *triout = 1;
                    break;
                case 2: //From triangle D to A, from a spin 0 to 2
                    *xout = xin;
                    *yout = yin + 1;
                    *triout = 0;
                    break;
                }
                break;

            case 1: //Handle the cases where we leave a triangle D from a spin 1
                switch (nnew) {
                case 0: //From triangle D to A, from a spin 1 to 0
                    *xout = xin;
                    *yout = yin;
                    *triout = 0;
                    break;
                case 2: //From triangle D to C, from a spin 1 to 2
                    *xout = xin - 1;
                    *yout = yin;
                    *triout = 2;
                    break;
                }
                break;

            case 2: //Handle the cases where we leave a triangle D from a spin 2
                switch (nnew) {
                case 0: //From triangle D to B, from a spin 2 to 0
                    *xout = xin;
                    *yout = yin;
                    *triout = 1;
                    break;
                case 1: //From triangle D to C, from a spin 2 to 1
                    *xout = xin;
                    *yout = yin;
                    *triout = 2;
                    break;
                }
                break;
            }
            break;
        } //end switch(triin)
    }   //end else(local=0)



    //We must now take in account the PBC.
    //Since we took a complicated cell made of 4 triangles (16 spins) the PBC are trivial.
    //We just need to reset the extrem value to the other extrem :
    if (L % 2 == 0) {
        if (*yout == L) {
            *yout = 0;
            *xout -= int(L / 2);
        }
        if (*yout == -1) {
            *yout = L - 1;
            *xout += int(L / 2);
        }
    }
    else if (L % 2 == 1) {
        if (*yout == L) {
            *yout = 0;
            *xout -= int(L / 2);
            switch (*triout) {
            case 0:
                *triout = 1; *xout -= 1;break;
            case 1:
                *triout = 0;            break;
            case 2:
                *triout = 3;            break;
            case 3:
                *triout = 2; *xout -= 1;break;
            }
        }
        if (*yout == -1) {
            *yout = L - 1;
            *xout += int(L / 2);
            switch (*triout) {
            case 0:
                *triout = 1;            break;
            case 1:
                *triout = 0; *xout += 1;break;
            case 2:
                *triout = 3; *xout += 1;break;
            case 3:
                *triout = 2;            break;
            }
        }

    }


    if (*xout >= L) { *xout -= L; }

    if (*xout <  0) { *xout += L ; }



}

int findSpinBetweenTwoOthers(const char tri1, const int n1, const char tri2)
{
    switch (tri1)
    {
    case 'A':
        switch (tri2)
        {
        case 'B':
            if (n1 == 1) { return(2); }
            if (n1 == 2) { return(1); }
            return(2);
            break;
        case 'C':
            if (n1 == 0) { return(2); }
            if (n1 == 1) { return(0); }
            break;
        case 'D':
            if (n1 == 0) { return(1); }
            if (n1 == 2) { return(0); }
            break;
        }
        break;

    case 'B':
        switch (tri2)
        {
        case 'A':
            if (n1 == 1) { return(2); }
            if (n1 == 2) { return(1); }
            return(2);
            break;
        case 'C':
            if (n1 == 0) { return(1); }
            if (n1 == 2) { return(0); }
            break;
        case 'D':
            if (n1 == 0) { return(2); }
            if (n1 == 1) { return(0); }
            break;
        }
        break;

    case 'C':
        switch (tri2)
        {
        case 'A':
            if (n1 == 0) { return(1); }
            if (n1 == 2) { return(0); }
            return(2);
            break;
        case 'B':
            if (n1 == 0) { return(2); }
            if (n1 == 1) { return(0); }
            break;
        case 'D':
            if (n1 == 1) { return(2); }
            if (n1 == 2) { return(1); }
            break;
        }
        break;

    case 'D':
        switch (tri2)
        {
        case 'A':
            if (n1 == 0) { return(2); }
            if (n1 == 1) { return(0); }
            return(2);
            break;
        case 'B':
            if (n1 == 0) { return(1); }
            if (n1 == 2) { return(0); }
            break;
        case 'C':
            if (n1 == 1) { return(2); }
            if (n1 == 2) { return(1); }
            break;
        }
        break;
    }

    std::cout << "There is an issue in findSpinBetweenTwoOthers : tri1 = " << tri1 << "   tri2 = " << tri2 << std::endl;
    return(-1);
}

bool isItFlipped(int x, int y, int tri, int b, int k, int tempLattice[L][L][4]) {
    //if the dimer didn't change, b didn't flip
    if (g_lattice[k][x][y][tri] == tempLattice[x][y][tri]) { return(false); }
    //if we are here, it means that the up triangle x, y tri have been flipped.
    //c flipped if and only if x, y, tri, b is a dimer of tempLattice or g_lattice.
    if (g_lattice[k][x][y][tri] == b || tempLattice[x][y][tri] == b) { return(true); }
    //if we are here, it means that the up triangle flipped from a to b (the two other possible position).
    return(false);
}