#include "Block.h"

Block::Block()
{
    id = -1;

    Nx = 0;
    Ny = 0;
}

Block::Block
(
    int id_,
    const Surface& surface_,
    int Nx_,
    int Ny_
)
{
    id = id_;

    surface = surface_;

    Nx = Nx_;
    Ny = Ny_;
}

void Block::generateMesh()
{
    mesh.generateCartesian
    (
        surface,
        Nx,
        Ny
    );
}