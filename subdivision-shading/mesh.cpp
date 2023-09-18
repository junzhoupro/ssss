#include "mesh.h"
#include "math.h"
#include "assert.h"
#include <sstream>

Mesh::Mesh() {
    qDebug() << "✓✓ Mesh constructor (Empty)";
    faceNormalsComputed = false;
    limitSupportSet = false;
}

Mesh::Mesh(OBJFile* loadedOBJFile) {
    faceNormalsComputed = false;
    level = 0;
    limitSupportSet = false;

    qDebug() << "✓✓ Mesh constructor (OBJ)";

    // Convert loaded OBJ file to HalfEdge mesh
    int numVertices, numHalfEdges, numFaces;
    int k, m, n;

    numVertices = loadedOBJFile->vertexCoords.size();
    numHalfEdges = 0;

    for (k=0; k<loadedOBJFile->faceValences.size(); k++) {
        numHalfEdges += loadedOBJFile->faceValences[k];
    }

    numFaces = loadedOBJFile->faceValences.size();

    // Note - resize() invokes the Vertex() constructor, reserve() does not.
    Vertices.reserve(numVertices);
    // If boundaries are present, reserve twice as much = worst case scenario
    HalfEdges.reserve(2*numHalfEdges);
    Faces.reserve(numFaces);

    // Add Vertices

    for (k=0; k<numVertices; k++) {
        // Coords (x,y,z), Out, Valence, Index
        Vertices.append(Vertex(loadedOBJFile->vertexCoords[k], 0, k));
        // Out and valence are unknown at this point.
    }

    qDebug() << "   # Vertices" << Vertices.capacity() << Vertices.size();

    int indexH = 0;
    int currentIndex = 0;

    // Initialize every entry of PotentialTwins with an empty QVector (using resize() )
    PotentialTwins.resize(loadedOBJFile->vertexCoords.size());

    // Add Faces and most of the HalfEdges

    for (m=0; m<numFaces; m++) {
        // Side, Val, Index
        Faces.append(Face(loadedOBJFile->faceValences[m], m));

        for (n=0; n<loadedOBJFile->faceValences[m]; n++) {
            // Target, Next, Prev, Twin, Poly, Index
            HalfEdges.append(HalfEdge(&Vertices[loadedOBJFile->faceCoordInd[currentIndex+n]],
                    nullptr,
                    nullptr,
                    nullptr,
                    &Faces[m],
                    indexH));

            // Next, Prev and Twin of the above HalfEdge have to be assigned later! Starting below...

            if (n > 0) {
                HalfEdges[indexH-1].next = &HalfEdges[indexH];
                HalfEdges[indexH].prev = &HalfEdges[indexH-1];

                // Append index of HalfEdge to list of OutgoingHalfEdges of its TailVertex.
                PotentialTwins[loadedOBJFile->faceCoordInd[currentIndex+n-1]].append(indexH);
            }
            indexH++;
        }

        // HalfEdges[indexH-1] is the most recent addition.
        Faces[m].side = &HalfEdges[indexH-1];

        HalfEdges[indexH-1].next = &HalfEdges[indexH-n];
        HalfEdges[indexH-n].prev = &HalfEdges[indexH-1];

        PotentialTwins[loadedOBJFile->faceCoordInd[currentIndex+n-1]].append(indexH-n);

        currentIndex += loadedOBJFile->faceValences[m];
    }

    qDebug() << "   # Faces" << Faces.capacity() << Faces.size();
    qDebug() << "   # HalfEdges" << HalfEdges.capacity() << HalfEdges.size();

    // Outs and Valences of vertices
    for (k=0; k<Vertices.size(); k++) {
        if (PotentialTwins[k].size() == 0) {
            qWarning() << " ! Isolated Vertex? PotentialTwins empty for Index" << k;
            dispVertInfo(&Vertices[k]);
            continue;
        }
        Vertices[k].out = &HalfEdges[PotentialTwins[k][0]];
        // Not the correct valence when on the boundary! Fixed below.
        Vertices[k].val = PotentialTwins[k].size();
    }

    setTwins(numHalfEdges, indexH);

    for (k=0; k<Vertices.size(); k++) {
        Vertex *v = &(Vertices[k]);
        v->boundary = vertOnBoundary(v);
        //qDebug() << "v" << k << v->val<< v->boundary << "BBBBBBBB" ;
    }

    PotentialTwins.clear();
    PotentialTwins.squeeze();

    makeEdgesSharp(loadedOBJFile->faceSharpness);
    //makeFacesSharp(2);

    qDebug() << "   # Updated HalfEdges" << HalfEdges.capacity() << HalfEdges.size();
}

// Assigns a sharpness value to an edge and updates the corresponding twin halfedge and connected vertices accordingly
void Mesh::makeEdgeSharp(HalfEdge* edge, unsigned int sharpness) {
    float sharp = std::max(float(sharpness), edge->sharpness);
    edge->sharpness = sharp;
    edge->twin->sharpness = sharp;
    edge->target->sharp = (sharp > 0 || edge->target->sharp);
    edge->twin->target->sharp = (sharp > 0 || edge->twin->target->sharp);
    edge->target->willSharp = (sharp > 1 || edge->target->willSharp);
    edge->twin->target->willSharp = (sharp > 1 || edge->twin->target->willSharp);
}

// Set the sharpness values of the mesh according to the data loaded from an obj or sharpness file
void Mesh::makeEdgesSharp(QVector<QVector<unsigned int>> sharpness) {
    HalfEdge* currentEdge, *firstEdge;
    Vertex* v;

    for (int i = 0; i < HalfEdges.size(); i++) {
        HalfEdges[i].sharpness = 0;
        HalfEdges[i].target->sharp = false;
        HalfEdges[i].target->willSharp = false;
    }

    for (int i = 0; i < std::min(sharpness.size(), Faces.size()); i++) {
        currentEdge = Faces[i].side;
        for (int j = 0; j < sharpness[i].size(); j++) {
            makeEdgeSharp(currentEdge, sharpness[i][j]);
            currentEdge = currentEdge->next;
        }
    }

    for (int i = 0; i < Vertices.size(); i++) {
        v = &Vertices[i];
        if (!v->sharp) {
            continue;
        }
        v->sharpMaskFaces = QVector<QVector<int>>(1);
        //v->sharpNormals = QVector<QVector3D>{QVector3D(1,0,0)};
        currentEdge = v->out;
        do {
            currentEdge = currentEdge->prev->twin;
        } while (currentEdge->sharpness == 0);
        firstEdge = currentEdge;
        do {
            v->sharpMaskFaces.last().append(currentEdge->polygon->index);
            currentEdge = currentEdge->prev->twin;
            if (firstEdge == currentEdge) {
                break;
            }
            if (currentEdge->sharpness > 0) {
                v->sharpMaskFaces.append(QVector<int>());
            }
        } while (true);
    }
}

void Mesh::setSharpnessFromFile(OBJFile sharpnessFile) {
    makeEdgesSharp(sharpnessFile.faceSharpness);
}

// Makes all edges of the first n faces sharp, just for debugging purposes
void Mesh::makeFacesSharp(int n) {

    for (int i = 0; i < HalfEdges.size(); i++) {
        HalfEdges[i].sharpness = 0;
        HalfEdges[i].target->sharp = false;
        HalfEdges[i].target->willSharp = false;
    }

    Vertex* v;
    HalfEdge* currentEdge, *firstEdge;
    for (int i = 0; i < n; i++) {
        currentEdge = Faces[i].side;
        firstEdge = currentEdge;

        do {
            currentEdge->sharpness = 2;
            currentEdge->twin->sharpness = 2;
            currentEdge->target->sharp = true;
            currentEdge->twin->target->sharp = true;
            currentEdge->target->willSharp = true;
            currentEdge->twin->target->willSharp = true;
            currentEdge = currentEdge->next;
        } while (currentEdge != firstEdge);
    }
    for (int i = 0; i < Vertices.size(); i++) {
        v = &Vertices[i];
        if (!v->sharp) {
            continue;
        }
        v->sharpMaskFaces = QVector<QVector<int>>(1);
        currentEdge = v->out;
        do {
            currentEdge = currentEdge->prev->twin;
        } while (currentEdge->sharpness == 0);
        firstEdge = currentEdge;
        do {
            v->sharpMaskFaces.last().append(currentEdge->polygon->index);
            currentEdge = currentEdge->prev->twin;
            if (firstEdge == currentEdge) {
                break;
            }
            if (currentEdge->sharpness > 0) {
                v->sharpMaskFaces.append(QVector<int>());
            }
        } while (true);
    }
}

// Checks if given vertex is vertex on boundary of the mesh
// And if so, it returns the corresponding edge
HalfEdge* Mesh::vertOnBoundary(Vertex* currentVertex) {

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

Mesh::~Mesh() {
    qDebug() << "✗✗ Mesh destructor";

    qDebug() << "   # Vertices:" << Vertices.size();
    qDebug() << "   # HalfEdges:" << HalfEdges.size();
    qDebug() << "   # Faces:" << Faces.size();

    Vertices.clear();
    Vertices.squeeze();
    HalfEdges.clear();
    HalfEdges.squeeze();
    Faces.clear();
    Faces.squeeze();
}

void Mesh::setTwins(int numHalfEdges, int indexH) {

    int m, n;
    int hTail, hHead, len;
    QSet<int> Twinless;

    // Assign Twins
    for (m=0; m<numHalfEdges; m++) {
        if (HalfEdges[m].twin == nullptr) {
            hTail = HalfEdges[m].prev->target->index;
            hHead = HalfEdges[m].target->index;
            len = HalfEdges[m].target->val;
            for (n=0; n<len; n++) {
                if (HalfEdges[PotentialTwins[hHead][n]].target->index == hTail) {
                    //qDebug() << "Found Twin!";
                    HalfEdges[m].twin = &HalfEdges[PotentialTwins[hHead][n]];
                    HalfEdges[PotentialTwins[hHead][n]].twin = &HalfEdges[m];
                    break;
                }
            }
            if (n == len) {
                // Twin not found...
                Twinless.insert(m);
            }
        }
    }

    if (Twinless.size() > 0) {
        qDebug() << " * There are" << Twinless.size() << "HalfEdges without Twin (i.e. the model contains boundaries)";
    }

    if (Twinless.size() > 0) {
        // The mesh is not closed

        //qDebug() << Twinless.values();

        HalfEdge* initialEdge;
        HalfEdge* currentEdge;
        int startBoundaryLoop;

        while (Twinless.size() > 0) {
            // Select a HalfEdge without Twin. The Twin that we will create is part of a boundary edge loop
            qDebug() << " → Processing new Boundary Edge Loop";

            initialEdge = &HalfEdges[*Twinless.begin()];
            Twinless.remove(initialEdge->index);

            // Target, Next, Prev, Twin, Poly, Index
            HalfEdges.append(HalfEdge( initialEdge->prev->target,
                                       nullptr,
                                       nullptr,
                                       initialEdge,
                                       nullptr,
                                       indexH ));
            startBoundaryLoop = indexH;
            // Twin of initialEdge should be assigned AFTER the central while loop!
            indexH++;

            // Use a sketch to properly understand these steps (assume counter-clockwise HalfEdges) :)
            currentEdge = initialEdge->prev;
            while (currentEdge->twin != nullptr) {
                currentEdge = currentEdge->twin->prev;
            }

            // Trace the current boundary loop
            while (currentEdge != initialEdge) {
                Twinless.remove(currentEdge->index);

                // Target, Next, Prev, Twin, Poly, Index
                HalfEdges.append(HalfEdge( currentEdge->prev->target,
                                           nullptr,
                                           &HalfEdges[indexH-1],
                                 currentEdge,
                                 nullptr,
                                 indexH ));
                HalfEdges[indexH-1].next = &HalfEdges[indexH];

                currentEdge->target->val += 1;
                currentEdge->twin = &HalfEdges[indexH];
                indexH++;

                currentEdge = currentEdge->prev;
                while (currentEdge->twin != nullptr) {
                    currentEdge = currentEdge->twin->prev;
                }
            }

            HalfEdges[startBoundaryLoop].prev = &HalfEdges[indexH-1];
            HalfEdges[indexH-1].next = &HalfEdges[startBoundaryLoop];

            initialEdge->target->val += 1;
            // Set Twin of initialEdge!
            initialEdge->twin = &HalfEdges[startBoundaryLoop];
        }

    }

}

// Compute normal of a face
void Mesh::setFaceNormal(Face* currentFace) {

    QVector3D faceNormal;
    HalfEdge *currentEdge;
    int k;
    double total;

    faceNormal = QVector3D();
    currentEdge = currentFace->side;

    // Sum cross products of edges and
    for (k=0; k<currentFace->val; k++) {
        faceNormal += 100.0*QVector3D::crossProduct(
                    currentEdge->next->target->origCoords - currentEdge->target->origCoords,
                    currentEdge->twin->target->origCoords - currentEdge->target->origCoords );
        currentEdge = currentEdge->next;
    }

    // Normalize normals sum
    faceNormal.normalize();

    if(!(faceNormal.length() > 0.99)) {
        total = abs(faceNormal.x()) + abs(faceNormal.y()) + abs(faceNormal.z());
        if(total != 0.0) {
            faceNormal *= 1.0 / total;

            faceNormal.normalize();

            if(!(faceNormal.length() > 0.99)) {
                qDebug() << "still too short faceNormal:" << faceNormal;
            }
        } else {
            qDebug() << "too short faceNormal:" << faceNormal;
        }
    }

    currentFace->normal = faceNormal;
}

QVector3D Mesh::normalizeNormal(QVector3D vertexNormal) {
    double total, fac;
    vertexNormal.normalize();
    //qDebug() << "vertexNormal" << vertexNormal << "vertexNormal.length()" << vertexNormal.length();

    if(!(vertexNormal.length() > 0.99)) {
        total = abs(vertexNormal.x()) + abs(vertexNormal.y()) + abs(vertexNormal.z());
        fac = 1.0 / total;

        if(total != 0.0) {
            vertexNormal *= fac;

            vertexNormal.normalize();

            if(!(vertexNormal.length() > 0.99)) {
                qDebug() << "still too short vertexNormal:" << vertexNormal << ",fac:" << fac;
            }

        }
        else {
            qDebug() << "too short vertexNormal:" << vertexNormal;
        }
    }

    return vertexNormal;
}

// Compute surface normal of a vertex by area-weighted facenormals
QVector3D Mesh::computeVertexSurfaceNormal(Vertex* currentVertex) {
    QVector3D vertexNormal;
    HalfEdge* currentEdge;
    QVector3D currentCoords, side1, side2;
    int k, normalIndex;
    double faceAngle, faceArea, total, fac;

    vertexNormal = QVector3D();
    currentEdge = currentVertex->out;
    QVector<QVector3D> vertexNormals(currentVertex->sharpMaskFaces.size());

    if (currentVertex->sharp) {//qDebug() << "SSSS";
        for (k=0; k<currentVertex->val; k++) {
            if (currentEdge->polygon) {
                for (int i = 0; i < currentVertex->sharpMaskFaces.size(); i++) {
                    if ((currentVertex->sharpMaskFaces[i].contains(currentEdge->polygon->parent))
                            || (currentVertex->sharpMaskFaces[i].contains(currentEdge->polygon->index)
                                && currentEdge->polygon->parent == -1)) {
                    //)){
                        normalIndex = i;
                        break;
                    }
                }//qDebug() << "normalIndex====" << normalIndex;
                currentCoords = currentVertex->origCoords;
                side1 = (currentEdge->target->origCoords - currentCoords);
                side2 = (currentEdge->prev->twin->target->origCoords - currentCoords);
                faceAngle = acos( fmax(-1.0, QVector3D::dotProduct(side1.normalized(), side2.normalized() ) ) );

                faceArea = (1000.0 * side1.length() * side2.length() * fmax(0.0, sin(faceAngle)));
                //qDebug() << "faceArea" << faceArea;
                //qDebug() << "currentEdge->polygon->normal" << currentEdge->polygon->normal;
                //qDebug() << "normalIndex" << normalIndex;
                vertexNormals[normalIndex] += faceArea * currentEdge->polygon->normal;
                //qDebug() << "vertexNormals[" << normalIndex << "]" << vertexNormals[normalIndex] ;
            }
            currentEdge = currentEdge->twin->next;
        }//for

        for (int i = 0; i < vertexNormals.size(); i++) {
            vertexNormals[i] = normalizeNormal(vertexNormals[i]);
        }
    }//sharp

    for (k=0; k<currentVertex->val; k++) {
        if (currentEdge->polygon) {
            currentCoords = currentVertex->origCoords;
            side1 = (currentEdge->target->origCoords - currentCoords);
            side2 = (currentEdge->prev->twin->target->origCoords - currentCoords);
            faceAngle = acos( fmax(-1.0, QVector3D::dotProduct(side1.normalized(), side2.normalized() ) ) );

            faceArea = (1000.0 * side1.length() * side2.length() * fmax(0.0, sin(faceAngle)));

            vertexNormal += faceArea * currentEdge->polygon->normal;
        }
        currentEdge = currentEdge->twin->next;
    }

    currentVertex->sharpSurfaceNormals = vertexNormals;
    vertexNormal.normalize();

//        if(!(vertexNormal.length() > 0.99)) {
//            total = abs(vertexNormal.x()) + abs(vertexNormal.y()) + abs(vertexNormal.z());
//            fac = 1.0 / total;
//            if(total != 0.0) {
//                vertexNormal *= fac;

//                vertexNormal.normalize();

//                if(!(vertexNormal.length() > 0.99)) {
//                    qDebug() << "still too short vertexNormal:" << vertexNormal << ",fac:" << fac;
//                }

//            } else {
//                qDebug() << "too short faceNormal:" << vertexNormal;
//            }}
    return vertexNormal;
}

// Compute surface normals for all faces
void Mesh::computeFaceNormals() {
    for (int k=0; k < Faces.size(); k++) {
        setFaceNormal(&Faces[k]);
    }
}

// Compute surface normals for all vertices
void Mesh::computeVertexSurfaceNormals() {
    Vertex *v;
    computeFaceNormals();

    for(int k = 0; k < Vertices.size(); k++) {
        v = &(Vertices[k]);
        v->surfaceNormal = computeVertexSurfaceNormal(v);
    }
    faceNormalsComputed = true;
}

// Set the normals of the vertices to the surface normals
void Mesh::setFaceNormals() {
    Vertex *v;

    if(!faceNormalsComputed) {
        computeVertexSurfaceNormals();
    }

    for(int k = 0; k < Vertices.size(); k++) {
        v = &(Vertices[k]);
        v->normal = v->surfaceNormal;
        v->sharpNormals = v->sharpSurfaceNormals;
        v->dotprodsSet = false;
    }
}

// Set the normals of the vertices to the subdivided normals
void Mesh::setSubdivideNormals() {
    Vertex *v;
    for(int k = 0; k < Vertices.size(); k++) {
        v = &(Vertices[k]);
        v->normal = v->sphereNormal;
        v->sharpNormals = v->sharpSphereNormals;
//        v->sharpNormals = v->sharpLinNormals;
//        qDebug() << "v.b" << v->boundary;
    }
}

// Set the normals of the vertices to the subdivided normals of a certain
// normal subdivision level (d = normalSubdivLevel - geomSubdivLevel)
void Mesh::setNormalsAtDepth(int d) {
    Vertex *v;
    for(int k = 0; k < Vertices.size(); k++) {
        v = &(Vertices[k]);
        v->normal = v->getNormalAtDepth(d);
        v->sharpNormals = v->getSharpNormalsAtDepth(d);
    }
}

// Set the normals of the vertices to the limit normals
void Mesh::setLimitNormals() {
    Vertex *v;
    for(int k = 0; k < Vertices.size(); k++) {
        v = &(Vertices[k]);
        v->normal = v->limitNormal;
        v->sharpNormals = v->sharpLimitNormals;
        if (v->sharpNormals.size() < v->sharpMaskVertices.size()) {
            v->sharpNormals.clear();
            for (int i = 0; i < v->sharpMaskVertices.size(); i++) {
                v->sharpNormals.append(v->limitNormal);
            }
        }
    }
}

// Blend the surface and the subdivided normals
void Mesh::blendNormals() {
    Vertex *v;
    for(int k = 0; k < Vertices.size(); k++) {
        v = &(Vertices[k]);
//        if(v->boundary)
//            qDebug() << "???????????" << v->boundary;
        //assert(v->limitBlendWeight <= 1.05);
        if (v->limitBlendWeight >= 1.05) {
            qDebug() << "Hmm, " << v->limitBlendWeight;
            v->limitBlendWeight = 1;
        }
        v->normal = (v->limitBlendWeight * v->sphereNormal + (1.0 - v->limitBlendWeight) * v->surfaceNormal).normalized();

        v->sharpNormals = QVector<QVector3D>(v->sharpSurfaceNormals.size());
        for (int i = 0; i < v->sharpNormals.size(); i++) {
            v->sharpNormals[i] = (v->limitBlendWeight * v->sharpSphereNormals[i] + (1.0 - v->limitBlendWeight) * v->sharpSurfaceNormals[i]).normalized();
        }
    }
}

// Set coordinates of the vertices to the original coordinates
void Mesh::setOrigCoords() {
    Vertex *v;
    for(int k = 0; k < Vertices.size(); k++) {
        v = &(Vertices[k]);
        v->coords = v->origCoords;
    }
}

// Set the coordinates of the vertices to the limit coordinates
void Mesh::setLimitCoords() {
    Vertex *v;
    for(int k = 0; k < Vertices.size(); k++) {
        v = &(Vertices[k]);
        v->coords = v->limitCoords;
    }
}

void Mesh::dispVertInfo(Vertex *dVert) {
    qDebug() << "Vertex at Index =" << dVert->index << "Coords =" << dVert->origCoords << "Out =" << dVert->out << "Val =" << dVert->val;
}

void Mesh::dispHalfEdgeInfo(HalfEdge *dHalfEdge) {
    qDebug() << "HalfEdge at Index =" << dHalfEdge->index << "Target =" << dHalfEdge->target << "Next =" << dHalfEdge->next << "Prev =" << dHalfEdge->prev << "Twin =" << dHalfEdge->twin << "Poly =" << dHalfEdge->polygon;
}

void Mesh::dispFaceInfo(Face *dFace){
    qDebug() << "Face at Index =" << dFace->index << "Side =" << dFace->side << "Val =" << dFace->val;
}

double Mesh::calcMinDotProduct() {
    double dotprod;
    double mindotprod = 1.0;
    for (int k = 0; k < Vertices.size(); k++) {
        Vertex *v = &(Vertices[k]);
        dotprod = v->dotprodLerpSlerp;
        mindotprod = fmin(mindotprod, dotprod);
    }
    qDebug() << "mindotprod: " << qSetRealNumberPrecision(20) << mindotprod;
    return mindotprod;
}

// Compute maximum difference between SLERP and LERP normals
// and return angle between corresponding normals
double Mesh::calcMaxSlerpAngle() {

    qDebug() << "acos(0.99999995) = " << qSetRealNumberPrecision( 10 ) << acos(0.99999995);

    return acos(calcMinDotProduct());
}

void Mesh::setSharpnessEdge(int edge, float sharpness) {
    HalfEdges[edge].sharpness = sharpness;
    HalfEdges[edge].twin->sharpness = sharpness;
}

void Mesh::substractEV(int meshType) {
    // meshType = 0: Loop
    // meshType = 1: C-C
    // Loop
    EVertices = QVector<Vertex*>();

    if (meshType == 0) {
        for (int k = 0; k < Vertices.size(); k++) {
            Vertex *v = &(Vertices[k]);
            // boundary && EV
            if(v->boundary && v->val != 3) {
                EVertices.append(v);
            }
            // not boundary && EV
            else if(v->val != 4) {
                EVertices.append(v);
            }

//            HalfEdge *currentEdge, *                                                                                                                                                                                                                                                                                                                                                                                ;
//            currentEdge = v->out;
//            for (int k=0; k<v->val; k++) {
//              // Vert in the star
//                EVertices.append(currentEdge->target);

//              faceEdge = currentEdge->next;

//              // Only the verts that are not in the star have to be updated.
//              // So not the vert itself, the star vert, or (after one round) the other star vert. Therefore, -3.
//              for (int l=0; l<faceEdge->polygon->val - 3; l++) {
//                  EVertices.append(faceEdge->target);
//                faceEdge = faceEdge->next;
//              }

//              currentEdge = currentEdge->twin->next;
//            }//
        }
    }
    // C-C
    else {
        for (int k = 0; k < Vertices.size(); k++) {
            Vertex *v = &(Vertices[k]);
            // boundary && EV
            if(v->boundary && v->val != 4) {
                EVertices.append(v);
            }
            // not boundary && EV
            else if(v->val != 6) {
                EVertices.append(v);
            }
        }
    }
//qDebug() << "EVertices" << EVertices;
}
