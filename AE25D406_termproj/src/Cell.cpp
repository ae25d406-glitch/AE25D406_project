#include "Cell.h"

#include <iostream>

Cell::Cell()
{
    id = -1;
}

Cell::Cell
(
    int id_,
    int n1,
    int n2,
    int n3,
    int n4
)
{
    id = id_;

    nodeIDs.push_back(n1);
    nodeIDs.push_back(n2);
    nodeIDs.push_back(n3);
    nodeIDs.push_back(n4);
}

void Cell::print() const
{
    std::cout << "Cell " << id << " : ";

    for (int node : nodeIDs)
    {
        std::cout << node << " ";
    }

    std::cout << std::endl;
}