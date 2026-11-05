#ifndef SURFACE_H
#define SURFACE_H

#include <vector>
#include "Line.h"

class Surface
{
public:

    int id;

    std::vector<Line*> boundaries;

    Surface();
    Surface(int id);

    void addBoundary(Line* line);

    void print() const;
};

#endif