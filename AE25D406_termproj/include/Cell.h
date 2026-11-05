#ifndef CELL_H
#define CELL_H

#include <vector>

class Cell
{
public:

    int id;

    std::vector<int> nodeIDs;

    Cell();

    Cell(
        int id_,
        int n1,
        int n2,
        int n3,
        int n4);

    void print() const;
};

#endif