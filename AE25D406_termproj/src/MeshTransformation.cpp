//------------------------------------------------------------//
// File: MeshTransformation.cpp
//
// Description:
//
// Applies geometric transformations to a mesh.
//
// Connectivity is NOT modified.
// Only node coordinates are updated.
//
//------------------------------------------------------------//

#include "MeshTransformation.h"

#include <cmath>

using namespace std;

//------------------------------------------------------------//
// Translation
//------------------------------------------------------------//

void MeshTransformation::translate
(
    Mesh& mesh,
    double dx,
    double dy
)
{
    for(auto& node : mesh.nodes)
    {
        node.x += dx;
        node.y += dy;
    }
}

//------------------------------------------------------------//
// Scaling
//------------------------------------------------------------//

void MeshTransformation::scale
(
    Mesh& mesh,
    double sx,
    double sy
)
{
    for(auto& node : mesh.nodes)
    {
        node.x *= sx;
        node.y *= sy;
    }
}

//------------------------------------------------------------//
// Rotation about Origin
//------------------------------------------------------------//

void MeshTransformation::rotate
(
    Mesh& mesh,
    double angleDeg
)
{
    double theta = angleDeg * M_PI / 180.0;

    double c = cos(theta);
    double s = sin(theta);

    for(auto& node : mesh.nodes)
    {
        double x = node.x;
        double y = node.y;

        node.x = c*x - s*y;
        node.y = s*x + c*y;
    }
}

//------------------------------------------------------------//
// Shearing
//------------------------------------------------------------//

void MeshTransformation::shear
(
    Mesh& mesh,
    double kx,
    double ky
)
{
    for(auto& node : mesh.nodes)
    {
        double x = node.x;
        double y = node.y;

        node.x = x + kx*y;
        node.y = y + ky*x;
    }
}

//------------------------------------------------------------//
// Sinusoidal Transformation
//------------------------------------------------------------//

void MeshTransformation::sinusoidal
(
    Mesh& mesh,
    double amplitude,
    double wavelength
)
{
    const double pi = acos(-1.0);

    for(auto& node : mesh.nodes)
    {
        node.y += amplitude *
                  sin(2.0*pi*node.x/wavelength);
    }
}