//------------------------------------------------------------//
// File: testGeometry.cpp
//
// Regression tests for
//
// Point
// Line
// Surface
//
//------------------------------------------------------------//

#include <iostream>

#include "Point.h"
#include "Line.h"
#include "Surface.h"

#include "Test.h"

using namespace std;

int main()
{
    cout << "=========================================\n";
    cout << " Running Geometry Tests\n";
    cout << "=========================================\n\n";

    Test test;

    //--------------------------------------------------------
    // Create Geometry
    //--------------------------------------------------------

    Point p0(0,0.0,0.0);
    Point p1(1,5.0,0.0);
    Point p2(2,5.0,2.0);
    Point p3(3,0.0,2.0);

    //--------------------------------------------------------

    Line l0(0,&p0,&p1);
    Line l1(1,&p1,&p2);
    Line l2(2,&p2,&p3);
    Line l3(3,&p3,&p0);

    //--------------------------------------------------------

    Surface rectangle(0);

    rectangle.addBoundary(&l0);
    rectangle.addBoundary(&l1);
    rectangle.addBoundary(&l2);
    rectangle.addBoundary(&l3);

    //--------------------------------------------------------
    // Point Tests
    //--------------------------------------------------------

    test.expectEqual(p0.x,0.0,1e-12,"Point 0 x");

    test.expectEqual(p0.y,0.0,1e-12,"Point 0 y");

    test.expectEqual(p2.x,5.0,1e-12,"Point 2 x");

    test.expectEqual(p2.y,2.0,1e-12,"Point 2 y");

    //--------------------------------------------------------
    // Line Tests
    //--------------------------------------------------------

    test.expectTrue
    (
        l0.start==&p0,
        "Line0 start point"
    );

    test.expectTrue
    (
        l0.end==&p1,
        "Line0 end point"
    );

    test.expectTrue
    (
        l3.end==&p0,
        "Line3 end point"
    );

    //--------------------------------------------------------
    // Surface Tests
    //--------------------------------------------------------

    test.expectEqual
    (
        rectangle.boundaries.size(),
        4,
        1e-12,
        "Surface has four boundaries"
    );

    test.expectEqual
    (
        rectangle.id,
        0,
        1e-12,
        "Surface ID"
    );

    //--------------------------------------------------------

    test.summary();

    return test.success() ? 0 : 1;
}