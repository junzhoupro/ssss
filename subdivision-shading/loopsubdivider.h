#ifndef LOOPSUBDIVIDER_H
#define LOOPSUBDIVIDER_H

#include "subdivider.h"

class LoopSubdivider : public Subdivider {
public:
    LoopSubdivider();

    // Implementation of Subdivider::subdivide
    void subdivide(Mesh *inputMesh, Mesh *subdivMesh);

    // Implementation of Subdivider::setLimitSupport
    void setLimitSupport(Mesh* mesh);

    // Implementation of Subdivider::setInitialBlendWeights
    void setInitialBlendWeights(Mesh *mesh, bool initZeros, short subdivApproach, int sharpEV, bool sharpEnd);

    // Implementation of Subdivider::setLinearBlendWeights
    void setLinearBlendWeights(Mesh *mesh, int sharpEV, bool sharpEnd);

private:
    // Adds a vertex point according to Loop subdiv. stencils
    Vertex vertexPoint(HalfEdge *firstEdge, short vval, int vindex);

    // Adds an edge point according to Loop subdiv. stencils
    Vertex edgePoint(HalfEdge *firstEdge, short vval, int vindex);

    // 3 functions used by subdivide() to create different types of vertices
    void createVertexPoints(Mesh *inputMesh, Mesh *subdivMesh);
    void createEdgePoints(Mesh *inputMesh, Mesh *subdivMesh, int numHalfEdges);
    void createFaces(Mesh *inputMesh, Mesh *subdivMesh, int numHalfEdges, int numFaces);

    // Determines the mask vertices and weights of a sharp vertex
    void setSharpMaskVertexPoint(Vertex* vertex, HalfEdge* firstEdge);
    void setSharpMaskEdgePoint(Vertex* vertex, HalfEdge* firstEdge);

    // Finds all sharp edges connected to a vertex
    QVector<HalfEdge*> getSharpEdges(Vertex* vertex);
};

#endif // LOOPSUBDIVIDER_H
