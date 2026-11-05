#ifndef MULTIBLOCKMESH_H
#define MULTIBLOCKMESH_H

#include <vector>

#include "Block.h"
#include "Mesh.h"

class MultiBlockMesh : public Mesh
{
public:

    std::vector<Block> blocks;

    MultiBlockMesh();

    void addBlock(const Block& block);

    void generate();

};

#endif