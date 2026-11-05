#ifndef MESHWRITER_H
#define MESHWRITER_H

#include <string>

#include "Mesh.h"

class MeshWriter
{

public:

    static void writeVTK
    (
        const Mesh& mesh,
        const std::string& filename
    );

};

#endif