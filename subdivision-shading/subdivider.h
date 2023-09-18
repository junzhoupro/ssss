#ifndef SUBDIVIDER_H
#define SUBDIVIDER_H

#include "mesh.h"

class Subdivider {

public:
    Subdivider();

    // function that subdivides inputMesh and writes result to subdivMesh
    virtual void subdivide(Mesh* inputMesh, Mesh* subdivMesh) = 0;

    // function that sets initial blend weights for given mesh and subdivision scheme
    virtual void setInitialBlendWeights(Mesh* mesh, bool initZeros, short subdivApproach, int sharpEV, bool sharpEnd) = 0;

    // function that sets initial (linear) blend weights for given mesh
    virtual void setLinearBlendWeights(Mesh* mesh, int sharpEV, bool sharpEnd) = 0;
    void setLinearBlendWeights(Mesh* mesh, short regularInner, short regularBoundary, int sharpEV, bool sharpEnd);

    // function that sets 'support' (control) vertices needed for calculating limit position/normal/blend weights
    virtual void setLimitSupport(Mesh *inputMesh) = 0;

    // function that linearly interpolates the given vectors using the given weights
    QVector3D linearInterpolation(QVector<QVector3D> inputNorms, QVector<double> weights);

    void splitHalfEdges(Mesh* inputMesh, Mesh* subdivMesh, unsigned int numHalfEdges, unsigned int numVertPts, unsigned int numFacePts);

    // function that calculates the normals based on the normals of the previous level and the control vertices
    void calculateNormals(Mesh *inputMesh, short slerpIterations, int blendWay, short subdivApproach, int sharpEV, int subdividedType);

    // function that calculates the coords based on the coords of the previous level and the control vertices
    void calculateCoords(Mesh *inputMesh);

    // function that calculates the (subdivision) blend weights based on the weights of the previous level and the control vertices
    void calculateBlendWeights(Mesh *inputMesh, short subdivApproach, int sharpEV, bool sharpEnd);

    // function that calculates the (linear) blend weights based on the weights of the previous level and the control vertices
    void calculateLinearBlendWeights(Mesh *inputMesh);

    // function that calculates the limit coords based on the vertices of the current mesh
    void calculateLimitCoords(Mesh *inputMesh);

    // function that calculates the limit normals based on the vertices of the current mesh
    void calculateLimitNormals(Mesh *inputMesh);

    // function that calculates the limit blend weights based on the vertices of the current mesh
    // p is the power to which the blend weight is raised
    void calculateLimitBlendWeights(Mesh *mesh, float p);

    void calculateNoBlendWeights(Mesh* inputMesh);


protected:
    short innerValency, boundaryValency;

    // Function that calculates the 3D spherical interpolation of the given normal vectors
    QVector3D calcSphereNormal(QVector3D linNormal, QVector<QVector3D> inputNorms, QVector<double> weights, short slerpIterations);

    // Checks if given vertex is vertex on boundary of the mesh
    // And if so, it returns the corresponding edge
    HalfEdge *vertOnBoundary(Vertex *currentVertex);

    // function returning a vertex point based on vertices of the previous level
    Vertex boundaryVertexPoint(HalfEdge *boundaryEdge, Vertex *parentVertex, short vval, int vindex);

    // function returning an edge point based on vertices of the previous level
    Vertex boundaryEdgePoint(HalfEdge *currentEdge, short vval, int vindex);

    // Set boundary limit support (cubic B-spline)
    void setBoundaryLimitSupport(HalfEdge *boundaryEdge, Vertex *v);
};

#endif // SUBDIVIDER_H




