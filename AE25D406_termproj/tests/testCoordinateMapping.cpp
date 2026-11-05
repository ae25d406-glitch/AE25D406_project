//------------------------------------------------------------//
// File: testCoordinateMapping.cpp
//
// Regression Tests
//
// Coordinate Mapping
//
//------------------------------------------------------------//

#include <iostream>

#include "Point.h"
#include "Line.h"
#include "Surface.h"

#include "Mesh.h"
#include "CoordinateMapping.h"

#include "Test.h"

using namespace std;

int main()
{
    cout << "=========================================\n";
    cout << " Running Coordinate Mapping Tests\n";
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
    // Original Mesh
    //--------------------------------------------------------

    Mesh original;

    original.generateCartesian
    (
        rectangle,
        20,
        10
    );

    //--------------------------------------------------------
    // Identity Mapping
    //--------------------------------------------------------

    Mesh identity = original;

    CoordinateMapping::identity(identity);

    test.expectEqual(identity.nodes[0].x,
                     original.nodes[0].x,
                     1e-12,
                     "Identity x");

    test.expectEqual(identity.nodes[0].y,
                     original.nodes[0].y,
                     1e-12,
                     "Identity y");

    //--------------------------------------------------------
    // Trapezoidal Mapping
    //--------------------------------------------------------

    Mesh trapezoid = original;

    CoordinateMapping::trapezoidal
    (
        trapezoid,
        1.5
    );

    test.expectTrue(
        trapezoid.nodes.back().x >
        original.nodes.back().x,
        "Trapezoidal expands top"
    );

    //--------------------------------------------------------
    // Sinusoidal Mapping
    //--------------------------------------------------------

    Mesh wave = original;

    CoordinateMapping::sinusoidal
    (
        wave,
        0.25
    );

    test.expectTrue(
        wave.nodes[100].y !=
        original.nodes[100].y,
        "Wave modifies y"
    );

    //--------------------------------------------------------
    // Polar Mapping
    //--------------------------------------------------------

    Mesh polar = original;

    CoordinateMapping::polar
    (
        polar,
        1.0,
        2.0
    );

    double r =
        sqrt
        (
            polar.nodes[0].x*polar.nodes[0].x +
            polar.nodes[0].y*polar.nodes[0].y
        );

    test.expectTrue
    (
        r>=1.0 && r<=2.0,
        "Polar radius"
    );

    //--------------------------------------------------------
    // Topology Preservation
    //--------------------------------------------------------

    test.expectEqual
    (
        polar.getNumberOfNodes(),
        original.getNumberOfNodes(),
        1e-12,
        "Node count preserved"
    );

    test.expectEqual
    (
        polar.getNumberOfCells(),
        original.getNumberOfCells(),
        1e-12,
        "Cell count preserved"
    );

    test.expectEqual
    (
        polar.cells[0].nodeIDs[0],
        original.cells[0].nodeIDs[0],
        1e-12,
        "Connectivity preserved"
    );

    //--------------------------------------------------------

    test.summary();

    return test.success() ? 0 : 1;
}