//------------------------------------------------------------//
// File: MeshWriter.cpp
//
// Description:
//
// Exports a structured mesh to Legacy VTK format.
//
//------------------------------------------------------------//

#include "MeshWriter.h"

#include <fstream>
#include <iostream>

using namespace std;

void MeshWriter::writeVTK
(
    const Mesh& mesh,
    const string& filename
)
{
    ofstream file(filename);

    if(!file.is_open())
    {
        cout << "Error opening " << filename << endl;
        return;
    }

    //--------------------------------------------------------
    // Header
    //--------------------------------------------------------

    file << "# vtk DataFile Version 3.0\n";
    file << "Cartesian Mesh\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n\n";

    //--------------------------------------------------------
    // Points
    //--------------------------------------------------------

    file << "POINTS "
         << mesh.getNumberOfNodes()
         << " float\n";

    for(const auto& node : mesh.nodes)
    {
        file << node.x
             << " "
             << node.y
             << " "
             << 0.0
             << "\n";
    }

    file << "\n";

    //--------------------------------------------------------
    // Cells
    //--------------------------------------------------------

    int nCells = mesh.getNumberOfCells();

    file << "CELLS "
         << nCells
         << " "
         << nCells*5
         << "\n";

    for(const auto& cell : mesh.cells)
    {
        file << "4 ";

        for(int node : cell.nodeIDs)
        {
            file << node << " ";
        }

        file << "\n";
    }

    file << "\n";

    //--------------------------------------------------------
    // Cell Types
    //--------------------------------------------------------

    file << "CELL_TYPES "
         << nCells
         << "\n";

    for(int i=0;i<nCells;i++)
    {
        file << "9\n";
    }

    file.close();

    cout << "\n";
    cout << "Mesh exported to "
         << filename
         << endl;
}