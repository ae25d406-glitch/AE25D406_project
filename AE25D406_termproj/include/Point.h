#ifndef POINT_H
#define POINT_H

#include <iostream>

class Point
{
public:
    int id;
    double x;
    double y;

    Point();
    Point(int id, double x, double y);

    void print() const;
};

#endif