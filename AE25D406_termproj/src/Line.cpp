#include "Line.h"
#include <iostream>

Line::Line()
{
    id = -1;
    start = nullptr;
    end = nullptr;
}

Line::Line(int id_, Point* p1, Point* p2)
{
    id = id_;
    start = p1;
    end = p2;
}

void Line::print() const
{
    std::cout
        << "Line "
        << id
        << " : Point "
        << start->id
        << " -> Point "
        << end->id
        << std::endl;
}