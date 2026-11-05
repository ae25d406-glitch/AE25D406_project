//------------------------------------------------------------//
// File: testCartesianMesh.cpp
//
// Regression Tests
//
// Cartesian Mesh Generation
//
//------------------------------------------------------------//

#include <iostream>

#include "Point.h"
#include "Line.h"
#include "Surface.h"

#include "Mesh.h"
#include "Test.h"

using namespace std;

int main()
{
    cout << "=========================================\n";
    cout << " Running Cartesian Mesh Tests\n";
    cout << "=========================================\n\n";

    Test test;

    //--------------------------------------------------------
    // Geometry
    //--------------------------------------------------------

    Point p0(0,0.0,0.0);
    Point p1(1,5.0,0.0);
    Point p2(2,5.0,2.0);
    Point p3(3,0.0,2.0);

    Line l0(0,&p0,&p1);
    Line l1(1,&p1,&p2);
    Line l2(2,&p2,&p3);
    Line l3(3,&p3,&p0);

    Surface rectangle(0);

    rectangle.addBoundary(&l0);
    rectangle.addBoundary(&l1);
    rectangle.addBoundary(&l2);
    rectangle.addBoundary(&l3);

    //--------------------------------------------------------
    // Generate Mesh
    //--------------------------------------------------------

    Mesh mesh;

    mesh.generateCartesian
    (
        rectangle,
        20,
        10
    );

    //--------------------------------------------------------
    // Mesh Statistics
    //--------------------------------------------------------

    test.expectEqual(mesh.getNumberOfNodes(),
                     231,
                     1e-12,
                     "Number of Nodes");

    test.expectEqual(mesh.getNumberOfCells(),
                     200,
                     1e-12,
                     "Number of Cells");

    test.expectEqual(mesh.Nx,
                     20,
                     1e-12,
                     "Nx");

    test.expectEqual(mesh.Ny,
                     10,
                     1e-12,
                     "Ny");

    test.expectEqual(mesh.dx,
                     0.25,
                     1e-12,
                     "dx");

    test.expectEqual(mesh.dy,
                     0.20,
                     1e-12,
                     "dy");

    //--------------------------------------------------------
    // Bounding Box
    //--------------------------------------------------------

    test.expectEqual(mesh.xmin,
                     0.0,
                     1e-12,
                     "xmin");

    test.expectEqual(mesh.xmax,
                     5.0,
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
    // First Node
    //--------------------------------------------------------

    test.expectEqual(mesh.nodes[0].x,
                     0.0,
                     1e-12,
                     "First Node x");

    test.expectEqual(mesh.nodes[0].y,
                     0.0,
                     1e-12,
                     "First Node y");

    //--------------------------------------------------------
    // Last Node
    //--------------------------------------------------------

    test.expectEqual(mesh.nodes.back().x,
                     5.0,
                     1e-12,
                     "Last Node x");

    test.expectEqual(mesh.nodes.back().y,
                     2.0,
                     1e-12,
                     "Last Node y");

    //--------------------------------------------------------
    // First Cell
    //--------------------------------------------------------

    test.expectEqual(mesh.cells[0].nodeIDs[0],
                     0,
                     1e-12,
                     "Cell0 Node0");

    test.expectEqual(mesh.cells[0].nodeIDs[1],
                     1,
                     1e-12,
                     "Cell0 Node1");

    test.expectEqual(mesh.cells[0].nodeIDs[2],
                     22,
                     1e-12,
                     "Cell0 Node2");

    test.expectEqual(mesh.cells[0].nodeIDs[3],
                     21,
                     1e-12,
                     "Cell0 Node3");

    //--------------------------------------------------------
    // Last Cell
    //--------------------------------------------------------

    Cell last = mesh.cells.back();

    test.expectTrue(last.id == 199,
                    "Last Cell ID");

    //--------------------------------------------------------

    test.summary();

    return test.success() ? 0 : 1;
}