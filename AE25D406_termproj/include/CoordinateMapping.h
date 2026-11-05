#ifndef COORDINATEMAPPING_H
#define COORDINATEMAPPING_H

#include "Mesh.h"
#include <string>

class CoordinateMapping
{
public:

    //------------------------------------------------------------
    // Identity Mapping
    //------------------------------------------------------------

    static void identity(Mesh& mesh);

    //------------------------------------------------------------
    // Trapezoidal Mapping
    //------------------------------------------------------------

    static void trapezoidal
    (
        Mesh& mesh,
        double topScale
    );

    //------------------------------------------------------------
    // Sinusoidal Mapping
    //------------------------------------------------------------

    static void sinusoidal
    (
        Mesh& mesh,
        double amplitude
    );

    //------------------------------------------------------------
    // Polar Mapping
    //------------------------------------------------------------

    static void polar
    (
        Mesh& mesh,
        double rInner,
        double rOuter
    );

    
    // Conformal O-grid around a body  (term project)
    //
    // Mesh must be a Cartesian mesh on [0,2*pi] x [0,etaMax]
    // (x = xi around the body, y = eta away from the wall).
    // zeta = zeta0 + R exp(eta) exp(i xi)      (log-polar)
    // z    = zeta + b^2/zeta                    (Joukowski)
    // shape: "circle" (diameter 1), "ellipse" or "joukowski"
    // (chord 1, centred at the origin)
   

    static void conformalBody
    (
        Mesh& mesh,//the rectangle mesh
        const std::string& shape,//"circle","ellipse" or"joukowski"
        double thickness//ellipse: thickness/chord (0.12); joukowski: ε (0.10)
    );

};


// Joukowski / ellipse / circle map used by conformalBody

class ConformalMap
{
public:
    std::string shape;//which shape
    double b, R, zeta0x, zeta0y, scale, shiftx;//map parameters

    ConformalMap(const std::string& shape, double thickness); //set the parameters

    void map(double xi, double eta, double& x, double& y) const;//(ξ,η) → (x,y)

    // eta at which the grid reaches radius Rfar (in chords)
    double etaMax(double Rfar) const; //η value at the far-field radius
};

#endif