#ifndef LINE_H
#define LINE_H

#include "Point.h"

class Line
{
public:
    int id;

    Point* start;
    Point* end;

    Line();
    Line(int id, Point* p1, Point* p2);

    void print() const;
};

#endif