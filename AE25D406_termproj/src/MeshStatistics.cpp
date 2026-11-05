//------------------------------------------------------------//
// File: MeshStatistics.cpp
//
// Description:
//
// Prints useful information about a mesh.
//
// Works for both:
//   1. Structured Cartesian Mesh
//   2. Multi Block Mesh
//
//------------------------------------------------------------//

#include "MeshStatistics.h"

#include <iostream>
#include <iomanip>

using namespace std;

void MeshStatistics::print(const Mesh& mesh)
{
    cout << "\n";
    cout << "=====================================\n";
    cout << "          Mesh Statistics\n";
    cout << "=====================================\n\n";

    //--------------------------------------------------------
    // Topology
    //--------------------------------------------------------

    cout << "Topology\n";
    cout << "--------\n";

    cout << left << setw(25)
         << "Number of Nodes"
         << " : "
         << mesh.getNumberOfNodes()
         << endl;

    cout << left << setw(25)
         << "Number of Cells"
         << " : "
         << mesh.getNumberOfCells()
         << endl;

    //--------------------------------------------------------
    // Geometry
    //--------------------------------------------------------

    cout << "\nGeometry\n";
    cout << "--------\n";

    cout << left << setw(25)
         << "xmin"
         << " : "
         << mesh.xmin
         << endl;

    cout << left << setw(25)
         << "xmax"
         << " : "
         << mesh.xmax
         << endl;

    cout << left << setw(25)
         << "ymin"
         << " : "
         << mesh.ymin
         << endl;

    cout << left << setw(25)
         << "ymax"
         << " : "
         << mesh.ymax
         << endl;

    //--------------------------------------------------------
    // Structured Mesh Information
    //--------------------------------------------------------

    if(mesh.Nx > 0 && mesh.Ny > 0)
    {
        cout << "\nStructured Mesh\n";
        cout << "---------------\n";

        cout << left << setw(25)
             << "Nx"
             << " : "
             << mesh.Nx
             << endl;

        cout << left << setw(25)
             << "Ny"
             << " : "
             << mesh.Ny
             << endl;

        cout << left << setw(25)
             << "dx"
             << " : "
             << mesh.dx
             << endl;

        cout << left << setw(25)
             << "dy"
             << " : "
             << mesh.dy
             << endl;
    }
    else
    {
        cout << "\nStructured Mesh\n";
        cout << "---------------\n";
        cout << "Resolution : Not Applicable (Multi Block Mesh)\n";
    }

    cout << "\n=====================================\n";
}