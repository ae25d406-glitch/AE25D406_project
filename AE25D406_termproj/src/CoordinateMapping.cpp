#include "CoordinateMapping.h"

#include <cmath>

using namespace std;

//------------------------------------------------------------//
// Identity Mapping
//------------------------------------------------------------//

void CoordinateMapping::identity(Mesh& mesh)
{
    // No change
}

//------------------------------------------------------------//
// Trapezoidal Mapping
//------------------------------------------------------------//

void CoordinateMapping::trapezoidal
(
    Mesh& mesh,
    double topScale
)
{
    double ymin = mesh.ymin;
    double ymax = mesh.ymax;

    for(auto& node : mesh.nodes)
    {
        double eta = (node.y-ymin)/(ymax-ymin);

        node.x *= (1.0 + eta*(topScale-1.0));
    }
}

//------------------------------------------------------------//
// Sinusoidal Mapping
//------------------------------------------------------------//

void CoordinateMapping::sinusoidal
(
    Mesh& mesh,
    double amplitude
)
{
    const double pi = acos(-1.0);

    double xmin = mesh.xmin;
    double xmax = mesh.xmax;

    double L = xmax-xmin;

    for(auto& node : mesh.nodes)
    {
        node.y += amplitude*
                  sin(2.0*pi*(node.x-xmin)/L);
    }
}

//------------------------------------------------------------//
// Polar Mapping
//------------------------------------------------------------//

void CoordinateMapping::polar
(
    Mesh& mesh,
    double rInner,
    double rOuter
)
{
    double xmin = mesh.xmin;
    double xmax = mesh.xmax;

    double ymin = mesh.ymin;
    double ymax = mesh.ymax;

    for(auto& node : mesh.nodes)
    {
        double xi =
            (node.x-xmin)/(xmax-xmin);

        double eta =
            (node.y-ymin)/(ymax-ymin);

        double r =
            rInner +
            eta*(rOuter-rInner);

        double theta =
            xi*2.0*acos(-1.0);

        node.x =
            r*cos(theta);

        node.y =
            r*sin(theta);
    }
}

//------------------------------------------------------------//
// Conformal map:log-polar + Joukowski  
//------------------------------------------------------------//

#include <complex>
#include <stdexcept>

ConformalMap::ConformalMap(const std::string& shape_, double thickness)
{
    shape = shape_;
    scale = 1.0;
    shiftx = 0.0;// default: no scaling
    zeta0x = 0.0;//no shift
    zeta0y = 0.0;//circle centred at the origin

    if(shape=="circle")              // diameter 1, identity map
    {
        b = 0.0;
        R = 0.5;
        return;
    }
    else if(shape=="ellipse")        // centred circle -> ellipse
    {
        b = 1.0;
        R = b*sqrt((1.0+thickness)/(1.0-thickness));
    }
    else if(shape=="joukowski")      // circle shifted left -> airfoil
    {
        double eps = thickness;// ε = 0.10
        b = 1.0;
        zeta0x = -eps*b;//center 
        R = (b-zeta0x) + 0.02*b;      // small rounding of the trailing edge
    }
    else
    {
        throw std::runtime_error("unknown shape "+shape);
    }

    // scale so that chord = 1 and mid-chord is at x = 0
    double xmin = 1e30, xmax = -1e30;//start with impossible values
    for(int q=0;q<20000;q++)//20000 points around the body
    {
        double x,y;
        //uniform oordinate in zeta 
        map(2.0*acos(-1.0)*q/20000.0,0.0,x,y);//point at ξ = 2π q/20000 on the wall(η = 0)
        xmin = min(xmin,x);//LE
        xmax = max(xmax,x);//TE
    }
    scale  = 1.0/(xmax-xmin);//chord becomes 1
    shiftx = -0.5*(xmin+xmax)*scale;//mid-chord moved to x = 0
}

void ConformalMap::map(double xi, double eta, double& x, double& y) const
{
    std::complex<double> zeta =
        std::complex<double>(zeta0x,zeta0y)
        + R*exp(eta)*std::complex<double>(cos(xi),sin(xi));

    std::complex<double> z = zeta;// circle: z = ζ
    if(b!=0.0) z += b*b/zeta; // ellipse/airfoil: z = ζ + b²/ζ  (Joukowski)

    x = scale*z.real() + shiftx;
    y = scale*z.imag();
}

double ConformalMap::etaMax(double Rfar) const
{
    return log(Rfar/(scale*R));
}

void CoordinateMapping::conformalBody
(
    Mesh& mesh,
    const std::string& shape,
    double thickness
)
{
    ConformalMap cm(shape,thickness);// compute b, R, ζ0, scale, shiftx

    for(auto& node : mesh.nodes)
    {
        double xi  = node.x;// the rectangle's x is ξ (around the body)
        double eta = node.y;// the rectangle's y is η (away from the wall)
        cm.map(xi,eta,node.x,node.y);//node.x, node.y now hold physical x, y
    }
}
