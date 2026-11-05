#include "Surface.h"
#include <iostream>

Surface::Surface()
{
    id = -1;
}

Surface::Surface(int id_)
{
    id = id_;
}

void Surface::addBoundary(Line* line)
{
    boundaries.push_back(line);
}

void Surface::print() const
{
    std::cout << "\nSurface " << id << std::endl;

    std::cout << "Boundary Lines :" << std::endl;

    for (auto line : boundaries)
    {
        std::cout
            << "   Line "
            << line->id
            << " ("
            << line->start->id
            << " -> "
            << line->end->id
            << ")"
            << std::endl;
    }
}