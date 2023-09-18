#ifndef VERTEX
#define VERTEX

#include <QVector3D>
#include <QDebug>

#include "face.h"

// Forward declaration
class HalfEdge;

class Vertex {
public:
    QVector3D coords, origCoords, limitCoords;
    QVector3D normal, surfaceNormal, linNormal, sphereNormal, limitNormal, leftNormal, rightNormal;
    double dotprodLerpSlerp, dotprodLimit;
    HalfEdge* out;
    short val;
    int index;
    bool dotprodsSet;
    bool boundary;
    float origBlendWeight;
    float limitBlendWeight;
    Vertex *parent;
    QVector<Vertex*> limitVertices;
    QVector<double> limitWeights;
    QVector<Vertex*> linblendParents;
    bool sharp, willSharp, sharpEnd;
    int type; //vertexpoint, edgepoint, facepoint

    QVector<QVector<Vertex*>> sharpMaskVertices;
    QVector<QVector<double>> sharpMaskWeights;
    QVector<QVector<int>> sharpMaskFaces;
    QVector<QVector<Face*>> sharpMaskSides;
    QVector<float> sharpBlendWeights;
    QVector<QVector3D> sharpNormals, sharpSurfaceNormals, sharpLinNormals, sharpSphereNormals, sharpLimitNormals, sharpleftNormals, sharpRightNormals;
    QVector<double> sharpDotprodLerpSlerps;
    QVector<QVector3D> inputNorms;

    QVector<Vertex*> getTwoRings() {
        QVector<Vertex*> res = QVector<Vertex*>();
        HalfEdge *faceEdge, *currentEdge;
        currentEdge = this->out;

//        // Set everything in the 1-ring nbh to 1 (not just the star, which would be the loop above...)
//        for (int k=0; k<val; k++) {
//          // Vert in the star
//            res.append(currentEdge->target);

//          faceEdge = currentEdge->next;

//          // Only the verts that are not in the star have to be updated.
//          // So not the vert itself, the star vert, or (after one round) the other star vert. Therefore, -3.
//          for (int l=0; l<faceEdge->polygon->val - 3; l++) {
//            faceEdge->target->origBlendWeight = 1;
//            faceEdge = faceEdge->next;
//          }

//          currentEdge = currentEdge->twin->next;
//        }
    }

    // Inline constructors
    Vertex() {
        setDefaults();
    }

    Vertex(QVector3D vcoords, short vval, int vindex) {
        setDefaults();
        origCoords = vcoords;
        coords = origCoords;
        val = vval;
        index = vindex;
        sharp = false;
    }

    Vertex(short vval, int vindex, bool boundVertex, QVector<Vertex*> vertices, QVector<double> weights, QVector<Face*> faces, int vtype, Vertex* vparent = nullptr, bool vsharp = false) {
        setDefaults();
        val = vval;
        index = vindex;
        maskVertices = vertices;
        maskWeights = weights;
        maskFaces = faces;
        parent = vparent;
        boundary = boundVertex;
        sharp = vsharp;
        type = vtype;
        inputNorms = QVector<QVector3D>(vertices.size());
        for (int i = 0; i < vertices.size(); i++) {
            if (vertices[i]->willSharp) {
                if (faces[i]->parent == -1) {
                    inputNorms[i] = vertices[i]->getSharpNormal(faces[i]->index);
                } else {
                    inputNorms[i] = vertices[i]->getSharpNormal(faces[i]->parent);
                }
            } else {
                inputNorms[i] = vertices[i]->normal.normalized();
            }
        }
    }

    void setSharpMasks(QVector<QVector<Vertex*>> vertices, QVector<QVector<double>> weights, QVector<QVector<int>> faces, QVector<QVector<Face*>> sides) {
        sharpMaskWeights = weights;
        sharpMaskVertices = vertices;
        sharpMaskFaces = faces;
        sharpMaskSides = sides;
    }

    QVector<Vertex*> getMaskVertices() {
        return maskVertices;
    }

    QVector<Face*> getMaskFaces() {
        return maskFaces;
    }

    QVector<double> getMaskWeights() {
        return maskWeights;
    }

    QVector3D getSharpNormal(int face) {
        if (!sharp || face < 0) {
            return normal.normalized();
        }

        for (int i = 0; i < sharpMaskFaces.size(); i++) {
//            qDebug() << "i" << i;
//            qDebug() << "sharpMaskFaces[i].contains(face)" << sharpMaskFaces[i] << face;
//            qDebug() << "sharpNormals.size()" << sharpNormals.size() << "\n";
            if (sharpMaskFaces[i].contains(face) && sharpNormals.size() > i) {
                return sharpNormals[i].normalized();
            }
        }

        qDebug() << "Face index error........." << face << sharpMaskFaces;
//        for (int i = 0; i < sharpMaskFaces.size(); i++) {
//            qDebug() << "i" << i;
//            qDebug() << "sharpMaskFaces[i].contains(face)" << sharpMaskFaces[i] << face;
//            qDebug() << "sharpNormals.size()" << sharpNormals.size() << "\n";
//            if (sharpMaskFaces[i].contains(face) && sharpNormals.size() > i) {
//                return sharpNormals[i].normalized();
//            }
//        }
        return normal.normalized();
    }

    float getSharpBlendWeight(int face) {
        if (!sharp || face < 0) {
            return origBlendWeight;
        }

        for (int i = 0; i < sharpMaskFaces.size(); i++) {
            if (sharpMaskFaces[i].contains(face) && sharpBlendWeights.size() > i) {
                return sharpBlendWeights[i];
            }
        }

        //qDebug() << "Face index error for blend weights..." << face << sharpBlendWeights;
        return origBlendWeight;
    }

    // get dth descendant normal of this vertex
    // 'descendants' structure used if normal subdiv. level higher than geometry subdiv. level
    QVector3D getNormalAtDepth(int d) {
        return children[d]->normal;
    }

    QVector<QVector3D> getSharpNormalsAtDepth(int d) {
        return children[d]->sharpNormals;
    }

    // calls addChild()
    // adds this vertex to its predecessors
    void addToParents() {
        parent->addChild(this);
    }

private:
    QVector<Vertex*> maskVertices;
    QVector<double> maskWeights;
    QVector<Face*> maskFaces;
    QVector<Vertex*> children;

    // sets default values of object
    // used by several constructors
    void setDefaults() {
        coords = QVector3D();
        origCoords = QVector3D();
        limitCoords = QVector3D();
        normal = QVector3D();
        limitNormal = QVector3D();
        linNormal = QVector3D();
        surfaceNormal = QVector3D();
        sphereNormal = QVector3D();
        leftNormal = QVector3D();
        rightNormal = QVector3D();
        sharpNormals = QVector<QVector3D>();
        dotprodLerpSlerp = -1.0;
        dotprodLimit = -1.0;
        out = nullptr;
        val = -1;
        index = -1;
        dotprodsSet = false;
        parent = nullptr;
        children.reserve(10);
        children.clear();
        children.append(nullptr);
        boundary = false;
        origBlendWeight = 0.0;
        limitBlendWeight = 0.0;
        sharp = false;
        willSharp = false;
        type = 0;
        sharpEnd = false;
    }

    // adds descendants recursively
    // if vertex has parent, then first add descendant to parent
    // after adding to parents, add to vertex itself
    // every 'generation' contains a single vertex
    // 'descendants' structure used if normal subdiv. level higher than geometry subdiv. level
    void addChild(Vertex *child) {
        if(parent != nullptr) {
            parent->addChild(child);
        }
        children.append(child);
    }

};

#endif // VERTEX
