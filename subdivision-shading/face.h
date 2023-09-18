#ifndef FACE
#define FACE

#include <QVector3D>

// Forward declaration
class HalfEdge;

class Face {

public:
    HalfEdge* side;
    unsigned short val;
    unsigned int index, parent;
    QVector3D normal, testSharpNormal;

    // Inline constructors

    Face() {
        side = nullptr;
        val = 0;
        index = 0;
        normal = QVector3D();
        parent = -1;
    }

    Face(unsigned short fval, unsigned int findex, unsigned int par = -1) {
        side = nullptr;
        val = fval;
        index = findex;
        normal = QVector3D();
        parent = par;
    }
};

#endif // FACE
