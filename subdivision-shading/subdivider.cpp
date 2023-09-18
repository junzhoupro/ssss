#include "subdivider.h"
#include "math.h"
#include "assert.h"

Subdivider::Subdivider() {
    qDebug() << "✓✓ Subdivider constructor (Empty)";
}

// Function that rotates a vector around a given axis
// http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.pdf
QVector3D rotateVectorAroundAxis(QVector3D point, QVector3D axis, float angle) {
    float cosval = cos(angle);
    float sinval = sin(angle);
    point.normalize();
    axis.normalize();
    float dotprod = QVector3D::dotProduct(point, axis);
    float onemincos = 1-cosval;
    float fac1 = dotprod * onemincos;
    float x,y,z,u,v,w;
    QVector3D rotated = QVector3D();
    x = point.x();
    y = point.y();
    z = point.z();
    u = axis.x();
    v = axis.y();
    w = axis.z();
    rotated.setX(u*fac1 + x*cosval + (-w*y + v*z) * sinval);
    rotated.setY(v*fac1 + y*cosval + (w*x - u*z) * sinval);
    rotated.setZ(w*fac1 + z*cosval + (-v*x + u*y) * sinval);
    return rotated.normalized();
}

// Function that given a collection of vectors and weights, computes the linear interpolation
// between the given vectors
QVector3D Subdivider::  linearInterpolation(QVector<QVector3D> inputNorms, QVector<double> weights) {
    int cnt = weights.size();
    QVector3D output = QVector3D(0.0, 0.0, 0.0);

    for(int i = 0; i < cnt; i++) {
        output += weights[i] * inputNorms[i];
    }

    return output;
}

// Function that calculates the 3D spherical interpolation of the given normal vectors
QVector3D Subdivider::calcSphereNormal(QVector3D linNormal, QVector<QVector3D> inputNorms, QVector<double> weights, short slerpIterations) {
    int cnt = inputNorms.size();
    QVector3D sphereNorm, newNorm;
    QVector3D linCombExpMap, inpNorm;
    double angle;
    int i, j;

    // step 1: as a starting normal the linear combination of input normal vectors
    sphereNorm = linNormal;

    for(i = 0; i < slerpIterations; i++) {
        // step 2 and 3: exponential map for each input norm. and linear combination of these maps
        linCombExpMap = QVector3D(0.0,0.0,0.0);
        for(j = 0; j < cnt; j++) {
            inpNorm = inputNorms[j];
            angle = acos(fmin(1.0, QVector3D::dotProduct(inpNorm, sphereNorm)));
            linCombExpMap += weights[j] * angle * (inpNorm - (QVector3D::dotProduct(inpNorm, sphereNorm) * sphereNorm)).normalized();
        }
        // step 4: rotate old normal along axis defined by cross product between old and expmap
        newNorm = rotateVectorAroundAxis(sphereNorm, QVector3D::crossProduct(sphereNorm, linCombExpMap), linCombExpMap.length());
        sphereNorm = newNorm.normalized();
    }

    return sphereNorm;
}

// Calculate the LERP and SLERP subdivided normals and their differences
void Subdivider::calculateNormals(Mesh* inputMesh, short slerpIterations, int blendWay, short subdivApproach, int sharpEV, int subdividedType) {
    QVector<Vertex*> supportVerts;
    QVector<QVector3D> inputNorms;
    QVector<double> weights;
    inputNorms.reserve(10);
    inputMesh->computeVertexSurfaceNormals();
    int wrongNormCnt = 0;

    // Improve efficiency

    // normal way
//    else {

    for(int i = 0; i < inputMesh->Vertices.size(); i++) {
        Vertex *v = &(inputMesh->Vertices[i]);
//        if(v->boundary)
//            qDebug() << "vb" << v->boundary;

        //not sharp
//        if (!v->sharp){

            inputNorms.clear();
                    supportVerts = v->getMaskVertices();
                    for(int j = 0; j < supportVerts.size(); j++) {
                        inputNorms.append(supportVerts[j]->normal);
                    }

        weights = v->getMaskWeights();
        v->linNormal = linearInterpolation(v->inputNorms, weights).normalized();
        wrongNormCnt += (v->linNormal.length() < 0.95);
        v->sphereNormal = calcSphereNormal(v->linNormal, v->inputNorms, weights, slerpIterations);
        v->dotprodLerpSlerp = QVector3D::dotProduct(v->linNormal, v->sphereNormal);
        v->dotprodsSet = true;
//        }

        //sharp
//        else
        if (v->sharp) {
            v->sharpLinNormals = QVector<QVector3D>(v->sharpMaskVertices.size());
            v->sharpSphereNormals = QVector<QVector3D>(v->sharpMaskVertices.size());
            v->sharpDotprodLerpSlerps = QVector<double>(v->sharpMaskVertices.size());

            for (int k = 0; k < v->sharpMaskVertices.size(); k++) {
                inputNorms.clear();
                for(int j = 0; j < v->sharpMaskVertices[k].size(); j++) {
                    int idx = v->sharpMaskSides[k][j]->parent;
                    if (idx == -1) {
                        idx = v->sharpMaskSides[k][j]->index;
                    }// if

                    inputNorms.append(v->sharpMaskVertices[k][j]->getSharpNormal(idx));
                }// for j
                weights = v->sharpMaskWeights[k];
                v->sharpLinNormals[k] = linearInterpolation(inputNorms, weights).normalized();
                wrongNormCnt += (v->linNormal.length() < 0.95);
                //Don't calculate subdivided normals when sharp
                v->sharpSphereNormals[k] = calcSphereNormal(v->sharpLinNormals[k], inputNorms, weights, slerpIterations);
//                v->sharpSphereNormals[k] = v->sharpLinNormals[k];
//                v->sharpSphereNormals[k] = QVector3D(0,0,0);
                v->sharpDotprodLerpSlerps[k] = QVector3D::dotProduct(v->sharpLinNormals[k], v->sharpSphereNormals[k]);
            }
        }

//        if(v->boundary)
//            qDebug() << "after vb" << v->boundary;
    }//for v
//    }
}

// Calculate the coordinates of the Mesh vertices given their control vertices
// and corresponding weights
void Subdivider::calculateCoords(Mesh* inputMesh) {
    QVector<Vertex*> supportVerts;
    QVector<QVector3D> inputCoords;
    QVector<double> weights;

    int boundaryev = 0;
    int ev3 = 0, ev5 =0, ev6 = 0;
    for(int i = 0; i < inputMesh->Vertices.size(); i++) {
        Vertex *v = &(inputMesh->Vertices[i]);

        inputCoords.clear();
        inputCoords.reserve(v->getMaskVertices().size());
        supportVerts = v->getMaskVertices();
        weights = v->getMaskWeights();
        for(int j = 0; j < supportVerts.size(); j++) {
            inputCoords.append(supportVerts[j]->origCoords);
        }
        v->origCoords = linearInterpolation(inputCoords, weights);
        if(v->boundary) {
            if(v->val != 3) {
                boundaryev ++;
            }
        }
        else {
            if(v->val == 3) {
                ev3 ++;
            }
            else if(v->val == 5) {
                ev5 ++;
            }
        }

//if(v->boundary)
//    qDebug() << "###" << v->boundary;

    }

    qDebug() << "ev3" << ev3 << "ev5" << ev5 << "ev6" << ev6;

}

// Function for Subdivision Blending
// Calculates blend weights by interpolating blend weights
// of control vertices using the subdivision stencils
void Subdivider::calculateBlendWeights(Mesh* inputMesh, short subdivApproach, int sharpEV, bool sharpEnd) {
    qDebug() << "calculateBlendWeights";

    //sharpEV = sharpoption::zero
    //enum sharpoption {regular, boundary. zero}

    //sharpEnd = dart
    QVector<Vertex*> supportVerts;
    QVector<double> weights;
    int countForEv = 0;
    int countForSharp = 0;

    for(int i = 0; i < inputMesh->Vertices.size(); i++) {
        Vertex *v = &(inputMesh->Vertices[i]);
        supportVerts = v->getMaskVertices();
        weights = v->getMaskWeights();
        v->origBlendWeight = 0.0;

        for(int j = 0; j < supportVerts.size(); j++) {
            v->origBlendWeight += weights[j] * supportVerts[j]->origBlendWeight;
        }

        if (v->origBlendWeight == 1) {
            countForEv ++;
        }

        if (v->sharp) {
            qDebug() << "sharp!";
            countForSharp ++;
            for (int k = 0; k < v->sharpMaskVertices.size(); k++) {
                weights = v->sharpMaskWeights[k];
                for(int j = 0; j < v->sharpMaskVertices[k].size(); j++) {
                    int idx = v->sharpMaskSides[k][j]->parent;
                    if (idx == -1) {
                        idx = v->sharpMaskSides[k][j]->index;
                    }
                    //v->sharpBlendWeights[k] += weights[j] * v->sharpMaskVertices[k][j]->getSharpBlendWeight(idx);
                }
            }
            qDebug() << "sharp!" << v->sharpMaskVertices.size() << v->sharpBlendWeights;
        }// if sharp
    }// for every v

    // Take care of EFs
    // Init 1r NBH
    if(subdivApproach == 2) {qDebug() << "subdivApproach == 2";
        setInitialBlendWeights(inputMesh, false, subdivApproach, sharpEV, sharpEnd);
    }

}

void Subdivider::calculateNoBlendWeights(Mesh* inputMesh) {
    for(int i = 0; i < inputMesh->Vertices.size(); i++) {
        Vertex *v = &(inputMesh->Vertices[i]);
        v->origBlendWeight = 0.0;
        v->limitBlendWeight = 0.0;
    }
}

// Function that performs blend weights by linearly interpolating
// blend weights of control vertices if not a vertex point
// (which already has a blend weight from parent vertex)
void Subdivider::calculateLinearBlendWeights(Mesh* inputMesh) {
    double val;

    for(int i = 0; i < inputMesh->Vertices.size(); i++) {
        Vertex *v = &(inputMesh->Vertices[i]);
        if(!v->parent) {
            val = v->linblendParents.size();
            v->origBlendWeight = 0.0;
            for(int j = 0; j < val; j++) {
                v->origBlendWeight += v->linblendParents[j]->origBlendWeight;
            }
            v->origBlendWeight /= val;
        } else {
            v->origBlendWeight = v->parent->origBlendWeight;
        }
        v->limitBlendWeight = v->origBlendWeight;

        if (v->limitBlendWeight >= 1.05) {
            v->limitBlendWeight = 1;
        }
    }
}

// Compute limit positions of vertices using limit stencil
void Subdivider::calculateLimitCoords(Mesh* mesh) {
    Vertex *v;
    QVector3D limitCoords;

    // if limit control vertices and weights not set yet
    if(!mesh->limitSupportSet) {
        qDebug() << "gonna set limit support";
        setLimitSupport(mesh);
    }

    for(int k = 0; k < mesh->Vertices.size(); k++) {
        v = mesh->Vertices[k].out->twin->target;
        limitCoords = QVector3D(0.0,0.0,0.0);

        for(short i = 0; i < v->limitVertices.size(); i++) {
            limitCoords += v->limitWeights[i] * v->limitVertices[i]->origCoords;
        }
        v->limitCoords = limitCoords;

    }
}

// Project blend weights (of Subdivision Blending) to limit blend weights using
// same limit stencil as used for coordinates
void Subdivider::calculateLimitBlendWeights(Mesh* mesh, float p) {
    qDebug() << "calculateLimitBlendWeights";
    Vertex *v;
    float limitweight;

    // if limit control vertices and weights not set yet
    if(!mesh->limitSupportSet) {
        setLimitSupport(mesh);
    }

    for(int k = 0; k < mesh->Vertices.size(); k++) {
        v = mesh->Vertices[k].out->twin->target;
        limitweight = 0.0;

        for(short i = 0; i < v->limitVertices.size(); i++) {
            limitweight += v->limitWeights[i] * v->limitVertices[i]->origBlendWeight;
        }
        v->limitBlendWeight = std::pow(limitweight, p);

        if (v->limitBlendWeight >= 1.05) {
            v->limitBlendWeight = 1;
        }
        //qDebug() << k << v->origBlendWeight << v->sharpBlendWeights << v->limitBlendWeight;
        //if(v->boundary)
          //  qDebug() << "@@" << v->boundary;
    }
}


// Compute limit normals using same limit stencils as for coordinates
void Subdivider::calculateLimitNormals(Mesh* mesh) {
    Vertex *v;

    // if limit control vertices and weights not set yet
    if(!mesh->limitSupportSet) {
        setLimitSupport(mesh);
    }

    for(int k = 0; k < mesh->Vertices.size(); k++) {
        v = mesh->Vertices[k].out->twin->target;
        v->limitNormal = QVector3D(0.0,0.0,0.0);
        for(short i = 0; i < v->limitVertices.size(); i++) {
            v->limitNormal += v->limitWeights[i] * v->limitVertices[i]->normal;
        }
        v->limitNormal.normalize();
        v->dotprodLimit = QVector3D::dotProduct(v->normal, v->limitNormal);
    }
}

void Subdivider::splitHalfEdges(Mesh* inputMesh, Mesh* subdivMesh, unsigned int numHalfEdges, unsigned int numVertPts, unsigned int numFacePts) {
    unsigned int k, m;
    unsigned int vIndex;
    HalfEdge* currentEdge;

    vIndex = numFacePts + numVertPts;

    for (k=0; k<numHalfEdges; k++) {
        currentEdge = &inputMesh->HalfEdges[k];
        m = currentEdge->twin->index;

        // Target, Next, Prev, Twin, Poly, Index
        subdivMesh->HalfEdges.append(HalfEdge(nullptr, nullptr, nullptr, nullptr, nullptr, 2*k));
        subdivMesh->HalfEdges.append(HalfEdge(nullptr, nullptr, nullptr, nullptr, nullptr, 2*k+1));

        if (k < m) {
            subdivMesh->HalfEdges[2*k].target = &subdivMesh->Vertices[ vIndex ];
            subdivMesh->HalfEdges[2*k+1].target = &subdivMesh->Vertices[ numFacePts + currentEdge->target->index ];
            vIndex++;
        }
        else {
            subdivMesh->HalfEdges[2*k].target = subdivMesh->HalfEdges[2*m].target;
            subdivMesh->HalfEdges[2*k+1].target = &subdivMesh->Vertices[ numFacePts + currentEdge->target->index ];

            // Assign Twins
            subdivMesh->HalfEdges[2*k].twin = &subdivMesh->HalfEdges[2*m+1];
            subdivMesh->HalfEdges[2*k+1].twin = &subdivMesh->HalfEdges[2*m];
            subdivMesh->HalfEdges[2*m].twin = &subdivMesh->HalfEdges[2*k+1];
            subdivMesh->HalfEdges[2*m+1].twin = &subdivMesh->HalfEdges[2*k];

            // Boundary edges are added last when importing a mesh, so their index will always be higher than their twins.
            if (!currentEdge->polygon) {
                subdivMesh->HalfEdges[2*k].next = &subdivMesh->HalfEdges[2*k+1];
                subdivMesh->HalfEdges[2*k+1].prev = &subdivMesh->HalfEdges[2*k];

                if (currentEdge > currentEdge->next) {
                    m = currentEdge->next->index;
                    subdivMesh->HalfEdges[2*k+1].next = &subdivMesh->HalfEdges[2*m];
                    subdivMesh->HalfEdges[2*m].prev = &subdivMesh->HalfEdges[2*k+1];
                }

                if (currentEdge > currentEdge->prev) {
                    m = currentEdge->prev->index;
                    subdivMesh->HalfEdges[2*k].prev = &subdivMesh->HalfEdges[2*m+1];
                    subdivMesh->HalfEdges[2*m+1].next = &subdivMesh->HalfEdges[2*k];
                }
            }
        }
    }

    for (k=0; k<numHalfEdges; k++) {
        currentEdge = &inputMesh->HalfEdges[k];
        subdivMesh->HalfEdges[2*k].sharpness = currentEdge->sharpness - 1;
        subdivMesh->HalfEdges[2*k].twin->sharpness = currentEdge->sharpness - 1;
    }

    // Note that Next, Prev and Poly are not yet assigned at this point.

}

// Checks if given vertex is vertex on boundary of the mesh
// And if so, it returns the corresponding edge
HalfEdge* Subdivider::vertOnBoundary(Vertex* currentVertex) {

    short n = currentVertex->val;
    int k;
    HalfEdge* currentEdge = currentVertex->out;

    // Loop over all incident edges and check if on boundary
    // If edge on boundary -> return that edge
    for (k=0; k<n; k++) {
        if (!currentEdge->polygon) {
            return currentEdge;
        }
        currentEdge = currentEdge->prev->twin;
    }

    return nullptr;
}

// Computes a vertex point on the boundary, given the incident edge that is on the boundary
// The same stencil is used for Catmull-Clark and Loop subdivision
Vertex Subdivider::boundaryVertexPoint(HalfEdge *boundaryEdge, Vertex *parentVertex, short vval, int vindex) {
    QVector<Vertex *> maskVertices;
    QVector<double> maskWeights;
    QVector<Face*> maskFaces;

    maskVertices.reserve(3);
    maskWeights.reserve(3);
    maskFaces.reserve(3);
    maskVertices.append(boundaryEdge->target);
    maskVertices.append(boundaryEdge->twin->target);
    maskVertices.append(boundaryEdge->prev->twin->target);
    maskWeights.append(1.0/8.0);
    maskWeights.append(6.0/8.0);
    maskWeights.append(1.0/8.0);
    maskFaces.append(boundaryEdge->polygon);
    maskFaces.append(boundaryEdge->twin->polygon);
    maskFaces.append(boundaryEdge->polygon);

    return Vertex(vval, vindex, parentVertex->boundary, maskVertices, maskWeights, maskFaces, 0, parentVertex);
}


// Computes a new edge point on the boundary, given the edge that is on the boundary
// The same stencil is used for Catmull-Clark and Loop subdivision
Vertex Subdivider::boundaryEdgePoint(HalfEdge *currentEdge, short vval, int vindex) {
    QVector<Vertex *> maskVertices;
    QVector<double> maskWeights;
    QVector<Face*> maskFaces;

    // Compute new edge point as average of two endpoints of the edge
    maskVertices.reserve(2);
    maskWeights.reserve(2);
    maskFaces.reserve(3);

    maskVertices.append(currentEdge->target);
    maskVertices.append(currentEdge->twin->target);

    maskWeights.append(1.0/2.0);
    maskWeights.append(1.0/2.0);

    maskFaces.append(currentEdge->polygon);
    maskFaces.append(currentEdge->twin->polygon);

    Vertex vert;
//    if(!sharp)
         vert = Vertex(vval, vindex, true, maskVertices, maskWeights, maskFaces, 1);
//    else
//        Vertex vert = Vertex(vval, vindex, false, maskVertices, maskWeights, maskFaces, 1);
    // set vertices from which blend weight should be averaged for linear blending
    vert.linblendParents = maskVertices;
    return vert;
}

// Set boundary limit support (cubic B-spline)
void Subdivider::setBoundaryLimitSupport(HalfEdge *boundaryEdge, Vertex *v) {
    v->limitVertices.reserve(3);
    v->limitVertices.append(v);
    v->limitVertices.append(boundaryEdge->target);
    v->limitVertices.append(boundaryEdge->prev->twin->target);

    v->limitWeights.reserve(3);
    v->limitWeights.append(4.0/6.0);
    v->limitWeights.append(1.0/6.0);
    v->limitWeights.append(1.0/6.0);
}

// Sets initial blend weights: 1 at irregular vertices and 0 at regular vertices
// Used by both Catmull-Clark and Loop subdividers
void Subdivider::setLinearBlendWeights(Mesh* mesh, short regularInner, short regularBoundary, int sharpEV, bool sharpEnd) {
    Vertex *v;
    for(int k = 0; k < mesh->Vertices.size(); k++) {
        v = mesh->Vertices[k].out->twin->target;

        if (v->boundary && !v->sharp) {
            v->origBlendWeight = (1.0 - (v->val == regularBoundary));
        } else {
            v->origBlendWeight = (1.0 - (v->val == regularInner));
        }

        v->limitBlendWeight = v->origBlendWeight;

        if (v->limitBlendWeight >= 1.05) {
            qDebug() << "Hmm, " << v->limitBlendWeight;
            v->limitBlendWeight = 1;
        }
    }
}
