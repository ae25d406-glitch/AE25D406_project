//------------------------------------------------------------//
// File: testMultiBlock.cpp
//
// Regression Tests
//
// Multi Block Mesh
//
//------------------------------------------------------------//

#include <iostream>

#include "Point.h"
#include "Line.h"
#include "Surface.h"

#include "Block.h"
#include "MultiBlockMesh.h"

#include "Test.h"

using namespace std;

int main()
{
    cout << "=========================================\n";
    cout << " Running Multi Block Tests\n";
    cout << "=========================================\n\n";

    Test test;

    //--------------------------------------------------------
    // Block 1
    //--------------------------------------------------------

    Point p0(0,0.0,0.0);
    Point p1(1,5.0,0.0);
    Point p2(2,5.0,2.0);
    Point p3(3,0.0,2.0);

    Line l0(0,&p0,&p1);
    Line l1(1,&p1,&p2);
    Line l2(2,&p2,&p3);
    Line l3(3,&p3,&p0);

    Surface s0(0);

    s0.addBoundary(&l0);
    s0.addBoundary(&l1);
    s0.addBoundary(&l2);
    s0.addBoundary(&l3);

    //--------------------------------------------------------
    // Block 2
    //--------------------------------------------------------

    Point p4(4,5.0,0.0);
    Point p5(5,10.0,0.0);
    Point p6(6,10.0,2.0);
    Point p7(7,5.0,2.0);

    Line l4(4,&p4,&p5);
    Line l5(5,&p5,&p6);
    Line l6(6,&p6,&p7);
    Line l7(7,&p7,&p4);

    Surface s1(1);

    s1.addBoundary(&l4);
    s1.addBoundary(&l5);
    s1.addBoundary(&l6);
    s1.addBoundary(&l7);

    //--------------------------------------------------------
    // Blocks
    //--------------------------------------------------------

    Block block0(0,s0,20,10);
    Block block1(1,s1,20,10);

    //--------------------------------------------------------
    // MultiBlock
    //--------------------------------------------------------

    MultiBlockMesh mesh;

    mesh.addBlock(block0);
    mesh.addBlock(block1);

    mesh.generate();

    //--------------------------------------------------------
    // Number of blocks
    //--------------------------------------------------------

    test.expectEqual(mesh.blocks.size(),
                     2,
                     1e-12,
                     "Number of Blocks");

    //--------------------------------------------------------
    // Global mesh size
    //--------------------------------------------------------

    test.expectEqual(mesh.getNumberOfNodes(),
                     462,
                     1e-12,
                     "Global Nodes");

    test.expectEqual(mesh.getNumberOfCells(),
                     400,
                     1e-12,
                     "Global Cells");

    //--------------------------------------------------------
    // Bounding Box
    //--------------------------------------------------------

    test.expectEqual(mesh.xmin,
                     0.0,
                     1e-12,
                     "xmin");

    test.expectEqual(mesh.xmax,
                     10.0,
                     1e-12,
                     "xmax");

    test.expectEqual(mesh.ymin,
                     0.0,
                     1e-12,
                     "ymin");

    test.expectEqual(mesh.ymax,
                     2.0,
                     1e-12,
                     "ymax");

    //--------------------------------------------------------
    // IDs
    //--------------------------------------------------------

    test.expectEqual(mesh.nodes.front().id,
                     0,
                     1e-12,
                     "First Node ID");

    test.expectEqual(mesh.nodes.back().id,
                     461,
                     1e-12,
                     "Last Node ID");

    test.expectEqual(mesh.cells.front().id,
                     0,
                     1e-12,
                     "First Cell ID");

    test.expectEqual(mesh.cells.back().id,
                     399,
                     1e-12,
                     "Last Cell ID");

    //--------------------------------------------------------
    // Node coordinates
    //--------------------------------------------------------

    test.expectEqual(mesh.nodes.front().x,
                     0.0,
                     1e-12,
                     "First Node x");

    test.expectEqual(mesh.nodes.front().y,
                     0.0,
                     1e-12,
                     "First Node y");

    test.expectEqual(mesh.nodes.back().x,
                     10.0,
                     1e-12,
                     "Last Node x");

    test.expectEqual(mesh.nodes.back().y,
                     2.0,
                     1e-12,
                     "Last Node y");

    //--------------------------------------------------------
    // Summary
    //--------------------------------------------------------

    test.summary();

    return test.success() ? 0 : 1;
}