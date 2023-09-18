#ifndef CATCLARKSUBDIVIDER_H
#define CATCLARKSUBDIVIDER_H

#include "subdivider.h"

class CatClarkSubdivider : public Subdivider {
public:
    CatClarkSubdivider();

    // Implementation of Subdivider::subdivide
    void subdivide(Mesh *inputMesh, Mesh *subdivMesh);

    // Implementation of Subdivider::setLimitSupport
    void setLimitSupport(Mesh* mesh);

    // Implementation of Subdivider::setInitialBlendWeights
    void setInitialBlendWeights(Mesh *mesh, bool initZeros, short subdivApproach, int sharpEV, bool sharpEnd);

    // Implementation of Subdivider::setLinearBlendWeights
    void setLinearBlendWeights(Mesh *mesh, int sharpEV, bool sharpEnd);
private:

    // Adds an edge point according to Catmull-Clark subdiv. stencils
    Vertex edgePoint(HalfEdge *firstEdge, Mesh *subdivMesh, short vval, int vindex);

    // Adds a vertex point according to Catmull-Clark subdiv. stencils
    Vertex vertexPoint(HalfEdge *firstEdge, Mesh *subdivMesh, short vval, int vindex);

    // Adds a face point according to Catmull-Clark subdiv. stencils
    Vertex facePoint(HalfEdge *firstEdge, short vval, int vindex);

    // Finds all sharp edges connected to a vertex
    QVector<HalfEdge*> getSharpEdges(Vertex* vertex);

    // Determines the mask vertices and weights of a sharp vertex
    void setSharpMaskVertexPoint(Vertex* vertex, HalfEdge* firstEdge, Mesh* mesh);
    void setSharpMaskEdgePoint(Vertex* vertex, HalfEdge* firstEdge, Mesh* mesh);
};

#endif // CATCLARKSUBDIVIDER_H
