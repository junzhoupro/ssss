#ifndef MESH_H
#define MESH_H

#include <QVector>

#include "vertex.h"
#include "face.h"
#include "halfedge.h"

#include "objfile.h"

class Mesh {

public:
    Mesh();
    Mesh(OBJFile *loadedOBJFile);
    ~Mesh();

    void setTwins(int numHalfEdges, int indexH);

    QVector<Vertex> Vertices;
    QVector<QVector3D> Normals;
    QVector<Face> Faces;
    QVector<HalfEdge> HalfEdges;

    QVector<QVector<int>> PotentialTwins;

    QVector<Vertex*> EVertices;

    int level;
    bool limitSupportSet;

    enum SubdivType {LOOP = 0, CATCLARK = 1, MODBUTTERFLY = 2};

    // For debugging
    void dispVertInfo(Vertex* dVert);
    void dispHalfEdgeInfo(HalfEdge* dHalfEdge);
    void dispFaceInfo(Face* dFace);

    void computeVertexSurfaceNormals();
    double calcMaxSlerpAngle();
    void setFaceNormals();
    void setSubdivideNormals();
    void setNormalsAtDepth(int d);
    void setOrigCoords();
    void setLimitCoords();
    void setLimitNormals();
    HalfEdge *vertOnBoundary(Vertex *currentVertex);
    void blendNormals();
    double calcMinDotProduct();
    QVector3D normalizeNormal(QVector3D vertexNormal);

    void setSharpnessEdge(int edge, float sharpness);
    void makeEdgeSharp(HalfEdge* edge, unsigned int sharpness);
    void makeEdgesSharp(QVector<QVector<unsigned int>> sharpness);
    void setSharpnessFromFile(OBJFile sharpnessFile);
    void makeFacesSharp(int n);
    void substractEV(int meshType);

private:
    bool faceNormalsComputed;
    void setFaceNormal(Face* currentFace);
    QVector3D computeVertexSurfaceNormal(Vertex* currentVertex);
    QVector3D computeVertexSurfaceNormal(Vertex* currentVertex, int sharpIndex);
    void computeFaceNormals();

};

#endif // MESH_H
