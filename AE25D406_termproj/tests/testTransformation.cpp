//------------------------------------------------------------//
// File: testTransformation.cpp
//
// Regression Tests
//
// Mesh Transformation
//
//------------------------------------------------------------//

#include <iostream>

#include "Point.h"
#include "Line.h"
#include "Surface.h"

#include "Mesh.h"
#include "MeshTransformation.h"

#include "Test.h"

using namespace std;

int main()
{
    cout << "=========================================\n";
    cout << " Running Mesh Transformation Tests\n";
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
    // Translation
    //--------------------------------------------------------

    Mesh translated = original;

    MeshTransformation::translate
    (
        translated,
        2.0,
        1.0
    );

    test.expectEqual(translated.nodes[0].x,
                     2.0,
                     1e-12,
                     "Translation x");

    test.expectEqual(translated.nodes[0].y,
                     1.0,
                     1e-12,
                     "Translation y");

    //--------------------------------------------------------
    // Scaling
    //--------------------------------------------------------

    Mesh scaled = original;

    MeshTransformation::scale
    (
        scaled,
        2.0,
        3.0
    );

    test.expectEqual(scaled.nodes.back().x,
                     10.0,
                     1e-12,
                     "Scaling x");

    test.expectEqual(scaled.nodes.back().y,
                     6.0,
                     1e-12,
                     "Scaling y");

    //--------------------------------------------------------
    // Rotation
    //--------------------------------------------------------

    Mesh rotated = original;

    MeshTransformation::rotate
    (
        rotated,
        90.0
    );

    test.expectEqual(rotated.nodes[0].x,
                     0.0,
                     1e-10,
                     "Rotation x");

    test.expectEqual(rotated.nodes[0].y,
                     0.0,
                     1e-10,
                     "Rotation y");

    //--------------------------------------------------------
    // Shear
    //--------------------------------------------------------

    Mesh sheared = original;

    MeshTransformation::shear
    (
        sheared,
        0.5,
        0.0
    );

    double expectedX =
        original.nodes.back().x +
        0.5*original.nodes.back().y;

    test.expectEqual(sheared.nodes.back().x,
                     expectedX,
                     1e-12,
                     "Shear x");

    //--------------------------------------------------------
    // Connectivity
    //--------------------------------------------------------

    test.expectEqual(translated.getNumberOfNodes(),
                     original.getNumberOfNodes(),
                     1e-12,
                     "Translation preserves nodes");

    test.expectEqual(translated.getNumberOfCells(),
                     original.getNumberOfCells(),
                     1e-12,
                     "Translation preserves cells");

    test.expectEqual(rotated.cells[0].nodeIDs[0],
                     original.cells[0].nodeIDs[0],
                     1e-12,
                     "Connectivity preserved");

    //--------------------------------------------------------

    test.summary();

    return test.success() ? 0 : 1;
}