#include "MultiBlockMesh.h"

#include <iostream>

using namespace std;

MultiBlockMesh::MultiBlockMesh()
{
}

void MultiBlockMesh::addBlock(const Block& block)
{
    blocks.push_back(block);
}

void MultiBlockMesh::generate()
{
    //----------------------------------------------------
    // Clear existing mesh
    //----------------------------------------------------

    nodes.clear();
    cells.clear();

    int globalNodeID = 0;
    int globalCellID = 0;

    //----------------------------------------------------
    // Loop over all blocks
    //----------------------------------------------------

    for(auto& block : blocks)
    {
        //--------------------------------------------
        // Generate block mesh
        //--------------------------------------------

        block.generateMesh();

        //--------------------------------------------
        // Store node offset
        //--------------------------------------------

        int nodeOffset = globalNodeID;

        //--------------------------------------------
        // Copy nodes
        //--------------------------------------------

        for(auto node : block.mesh.nodes)
        {
            node.id = globalNodeID++;

            nodes.push_back(node);
        }

        //--------------------------------------------
        // Copy cells
        //--------------------------------------------

        for(auto cell : block.mesh.cells)
        {
            cell.id = globalCellID++;

            for(auto& nodeID : cell.nodeIDs)
            {
                nodeID += nodeOffset;
            }

            cells.push_back(cell);
        }
    }
    //----------------------------------------------------
    // Compute Global Bounding Box
    //----------------------------------------------------

    computeBoundingBox();

    cout << "\n";
    cout << "=====================================\n";
    cout << "Multi Block Mesh Generated\n";
    cout << "=====================================\n";

    cout << "Blocks : "
         << blocks.size()
         << endl;

    cout << "Nodes : "
         << nodes.size()
         << endl;

    cout << "Cells : "
         << cells.size()
         << endl;
}