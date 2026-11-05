#ifndef OGRID_H
#define OGRID_H

// Finite-volume geometry of the conformal O-grid 
// Built from a Mesh produced by CoordinateMapping::conformalBody.
// Cell (I,J):  I = 0..Ni-1 around the body (periodic),
//              J = 0..Nj-1 away from the wall (J=0 touches the body).
// Flat index  k = I*Nj + J.
// Node (i,j) of the Mesh is nodes[j*(Ni+1)+i]

#include <vector>
#include <string>
#include "Mesh.h"

class OGrid
{
public:
    int Ni = 0, Nj = 0;

    std::vector<double> xc, yc, dV; // cell centres and cell volumes ΔV (area × unit depth)
    std::vector<double> Sex, Sey, gE;//(face vector) east face of cell k (towards I+1)
    std::vector<double> Snx, Sny, gN;// north face of cell k (towards J+1);
                                        // J=Nj-1: outer boundary, outward, g to the face
    std::vector<double> Swx, Swy, gW;   // wall face of (I,0), points into the body
    double bodyArea = 0;

    void build(const Mesh& mesh);//compute everything from the mapped mesh

    int n() const { return Ni*Nj; }//number of cells
    int id(int I,int J) const { return I*Nj+J; }//(I,J) → k
    int ip(int I) const { return (I+1==Ni) ? 0 : I+1; }//east neighbour column, wraps Ni−1 → 0
    int im(int I) const { return (I==0) ? Ni-1 : I-1; }// west neighbour column, wraps 0 → Ni−1
};

// Cartesian mesh on [0,2pi] x [0,etaMax(Rfar)] mapped by CoordinateMapping::conformalBody.
// xiStretch-moves grid lines away from the tips:  xi = s + a sin(2s)  (0 <= a < 0.5)
Mesh makeConformalMesh(const std::string& shape, double thickness, int Ni, int Nj, double Rfar,
                       double xiStretch = 0.0);

#endif
