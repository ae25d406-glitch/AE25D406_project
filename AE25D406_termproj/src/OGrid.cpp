#include "OGrid.h"
#include <cmath>

using namespace std;

static void orient(double& sx,double& sy,double dx,double dy)
{
    if(sx*dx+sy*dy<0.0){sx=-sx;sy=-sy;}//chek direction of face vector pointing out or in 
}

// g = |S| / (normal distance) = |S|^2 / (S . d)   (orthogonal grid)
static double gFace(double sx,double sy,double dx,double dy)
{
    return (sx*sx+sy*sy)/(sx*dx+sy*dy);
}

void OGrid::build(const Mesh& mesh)
{
    Ni = mesh.Nx;// cells around
    Nj = mesh.Ny; //cells outward 

    int N = Ni*Nj;
    //set side to N fills with zero 
    xc.assign(N,0); yc.assign(N,0); dV.assign(N,0);// one per cell
    Sex.assign(N,0); Sey.assign(N,0); Snx.assign(N,0); Sny.assign(N,0);// two faces per cell
    gE.assign(N,0); gN.assign(N,0);
    Swx.assign(Ni,0); Swy.assign(Ni,0); gW.assign(Ni,0);// one wall face per column

    auto X=[&](int i,int j){ return mesh.nodes[j*(Ni+1)+i].x; };// x of node (i,j)
    auto Y=[&](int i,int j){ return mesh.nodes[j*(Ni+1)+i].y; };// y of node (i,j)

    // centroids and areas (shoelace)
    for(int I=0;I<Ni;I++)
        for(int J=0;J<Nj;J++)
        {
            double px[4]={X(I,J),X(I+1,J),X(I+1,J+1),X(I,J+1)};// 4 corner in order
            double py[4]={Y(I,J),Y(I+1,J),Y(I+1,J+1),Y(I,J+1)};
            double a=0,cx=0,cy=0;
            for(int q=0;q<4;q++)//each edge
            {
                double cr=px[q]*py[(q+1)%4]-px[(q+1)%4]*py[q];// cross product x_q y_{q+1} − x_{q+1} y_q
                a+=cr;
                cx+=(px[q]+px[(q+1)%4])*cr;
                cy+=(py[q]+py[(q+1)%4])*cr;
            }
            int k=id(I,J);
            dV[k]=0.5*fabs(a);// area
            xc[k]=cx/(3.0*a);// centroid x
            yc[k]=cy/(3.0*a);// centroid y
        }

    // face area vectors
    for(int I=0;I<Ni;I++)
        for(int J=0;J<Nj;J++)
        {
            int k=id(I,J), kE=id(ip(I),J);// this cell and its east neighbour (wrapped)

            // east face: node (I+1,J) - (I+1,J+1)
            double sx= Y(I+1,J+1)-Y(I+1,J);//S =(Δy,−Δx )
            double sy=-(X(I+1,J+1)-X(I+1,J));
            orient(sx,sy,xc[kE]-xc[k],yc[kE]-yc[k]);//(center to center)point from P towards E
            Sex[k]=sx; Sey[k]=sy;// store S_e
            gE[k]=gFace(sx,sy,xc[kE]-xc[k],yc[kE]-yc[k]);//g_e = |S|²/(S·d_PE)

            // north face: node (I,J+1) -(I+1,J+1)
            sx= Y(I+1,J+1)-Y(I,J+1);// edge (I,J+1) -(I+1,J+1)
            sy=-(X(I+1,J+1)-X(I,J+1));
            if(J<Nj-1)// inner face
            {
                int kN=id(I,J+1);//north neighbour = k+1
                orient(sx,sy,xc[kN]-xc[k],yc[kN]-yc[k]);//point from P to N
                gN[k]=gFace(sx,sy,xc[kN]-xc[k],yc[kN]-yc[k]);// g_n = A/δx_PN
            }
            else
            {   //top row
                double fx=0.5*(X(I,J+1)+X(I+1,J+1)), fy=0.5*(Y(I,J+1)+Y(I+1,J+1)); // face midpoint
                orient(sx,sy,fx-xc[k],fy-yc[k]);//point outward (from P to the face)
                gN[k]=gFace(sx,sy,fx-xc[k],fy-yc[k]);//g = A/(distance P to face)
            }
            Snx[k]=sx; Sny[k]=sy;// store S_n

            // wall face: node (I,0)to (I+1,0), pointing into the body
            if(J==0)// only the cells touching the body
            {
                sx= Y(I+1,0)-Y(I,0);
                sy=-(X(I+1,0)-X(I,0));
                double fx=0.5*(X(I,0)+X(I+1,0)), fy=0.5*(Y(I,0)+Y(I+1,0));
                orient(sx,sy,fx-xc[k],fy-yc[k]);
                Swx[I]=sx; Swy[I]=sy;
                gW[I]=gFace(sx,sy,fx-xc[k],fy-yc[k]);
            }
        }

    // body area from the wall polygon
    double a=0;
    for(int I=0;I<Ni;I++) a+=X(I,0)*Y(I+1,0)-X(I+1,0)*Y(I,0);// shoelace over the wall nodes
    bodyArea=0.5*fabs(a);
}

#include "Point.h"
#include "Line.h"
#include "Surface.h"
#include "CoordinateMapping.h"

Mesh makeConformalMesh(const std::string& shape,double thickness,int Ni,int Nj,double Rfar,double xiStretch)
{
    ConformalMap cm(shape,thickness);// b, R, ζ0, scale for this shape
    double etaMax=cm.etaMax(Rfar);//ηmax = ln(Rfar/(scale·R))
    double twoPi=2.0*acos(-1.0);// 2π

    // computational rectangle (aquaMesh geometry classes)
    Point p0(0,0.0,0.0),p1(1,twoPi,0.0),p2(2,twoPi,etaMax),p3(3,0.0,etaMax);
    Line l0(0,&p0,&p1),l1(1,&p1,&p2),l2(2,&p2,&p3),l3(3,&p3,&p0);
    Surface s(0);
    s.addBoundary(&l0); s.addBoundary(&l1); s.addBoundary(&l2); s.addBoundary(&l3);

    Mesh mesh;
    mesh.generateCartesian(s,Ni,Nj); // uniform Ni * Nj grid on the rectangle
    // non-uniform spacing around the body: xi = s + a sin(2s)
    // (larger spacing at the tips xi = 0, pi; a < 0.5 keeps xi monotonic)
    for(auto& nd : mesh.nodes) nd.x += xiStretch*sin(2.0*nd.x);
    //in physical domain 
    CoordinateMapping::conformalBody(mesh,shape,thickness);
    return mesh;
}
