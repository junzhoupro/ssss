#include "catclarksubdivider.h"
#include "assert.h"

CatClarkSubdivider::CatClarkSubdivider() {
    qDebug() << "✓✓ CatClarkSubdivider constructor (Empty)";
}

// Implementation of Subdivider::subdivide() using the Catmull-Clark
// subdivision scheme
void CatClarkSubdivider::subdivide(Mesh* inputMesh, Mesh* subdivMesh) {
    int numVerts, numHalfEdges, numFaces, sumFaceValences;
    int k, m, s, t;
    int vIndex, hIndex, fIndex;
    short n;
    HalfEdge* currentEdge;

    qDebug() << ":: Creating new Catmull-Clark mesh";

    numVerts = inputMesh->Vertices.size();
    numHalfEdges = inputMesh->HalfEdges.size();
    numFaces = inputMesh->Faces.size();

    // Reserve memory
    subdivMesh->Vertices.reserve(numFaces + numVerts + numHalfEdges/2);

    sumFaceValences = 0;
    for (k=0; k<numFaces; k++) {
        sumFaceValences += inputMesh->Faces[k].val;
    }

    subdivMesh->HalfEdges.reserve(2*numHalfEdges + 2*sumFaceValences);
    subdivMesh->Faces.reserve(sumFaceValences);

    // Create face points
    for (k=0; k<numFaces; k++) {
        n = inputMesh->Faces[k].val;
        // Coords (x,y,z), Out, Valence, Index
        subdivMesh->Vertices.append( facePoint(inputMesh->Faces[k].side, n, k) );
    }

    qDebug() << " * Created face points";

    vIndex = numFaces;

    // Create vertex points
    for (k=0; k<numVerts; k++) {
        n = inputMesh->Vertices[k].val;
        // Coords (x,y,z), Out, Valence, Index
        subdivMesh->Vertices.append( vertexPoint(inputMesh->Vertices[k].out, subdivMesh, n, vIndex/*, inputMesh->Vertices[k].sharpness-1*/) );
        subdivMesh->Vertices[vIndex].addToParents();
        vIndex++;
    }

    qDebug() << " * Created vertex points";

    // Create edge points
    for (k=0; k<numHalfEdges; k++) {
        currentEdge = &inputMesh->HalfEdges[k];

        if (k < currentEdge->twin->index) {
            m = (!currentEdge->polygon || !currentEdge->twin->polygon) ? 3 : 4;
            // Coords (x,y,z), Out, Valence, Index

            Vertex v = edgePoint(currentEdge, subdivMesh, m, vIndex);
//            subdivMesh->Vertices.append( edgePoint(currentEdge, subdivMesh, m, vIndex) );
            subdivMesh->Vertices.append(v);
//            qDebug() << "new v" << v.boundary;
            vIndex++;
        }
    }

    qDebug() << " * Created edge points";
//    for(int i = 0; i < subdivMesh->Vertices.size(); i ++){
//        if(subdivMesh->Vertices[i].boundary)
//            qDebug() << "--wrong!!!!!";
//    }

    // Split halfedges
    splitHalfEdges(inputMesh, subdivMesh, numHalfEdges, numVerts, numFaces);

    qDebug() << " * Split halfedges";

    hIndex = 2*numHalfEdges;
    fIndex = 0;

    // Create faces and remaining halfedges
    for (k=0; k<numFaces; k++) {
        currentEdge = inputMesh->Faces[k].side;
        n = inputMesh->Faces[k].val;

        for (m=0; m<n; m++) {

            s = currentEdge->prev->index;
            t = currentEdge->index;

            // Side, Val, Index
            subdivMesh->Faces.append(Face(4, fIndex, inputMesh->Faces[k].index));

            subdivMesh->Faces[fIndex].side = &subdivMesh->HalfEdges[ 2*t ];

            // Target, Next, Prev, Twin, Poly, Index
            subdivMesh->HalfEdges.append(HalfEdge( &subdivMesh->Vertices[k],
                                                   nullptr,
                                                   &subdivMesh->HalfEdges[ 2*t ],
                                         nullptr,
                                         &subdivMesh->Faces[fIndex],
                                         hIndex ));

            subdivMesh->HalfEdges.append(HalfEdge( nullptr,
                                                   &subdivMesh->HalfEdges[2*s+1],
                                         &subdivMesh->HalfEdges[hIndex],
                                         nullptr,
                                         &subdivMesh->Faces[fIndex],
                                         hIndex+1 ));

            subdivMesh->HalfEdges[hIndex].next = &subdivMesh->HalfEdges[hIndex+1];
            subdivMesh->HalfEdges[hIndex+1].target = subdivMesh->HalfEdges[hIndex+1].next->twin->target;

            subdivMesh->HalfEdges[2*s+1].next = &subdivMesh->HalfEdges[2*t];
            subdivMesh->HalfEdges[2*s+1].prev = &subdivMesh->HalfEdges[hIndex+1];
            subdivMesh->HalfEdges[2*s+1].polygon = &subdivMesh->Faces[fIndex];

            subdivMesh->HalfEdges[2*t].next = &subdivMesh->HalfEdges[hIndex];
            subdivMesh->HalfEdges[2*t].prev = &subdivMesh->HalfEdges[2*s+1];
            subdivMesh->HalfEdges[2*t].polygon = &subdivMesh->Faces[fIndex];

            if (m > 0) {
                // Twins
                subdivMesh->HalfEdges[hIndex+1].twin = &subdivMesh->HalfEdges[hIndex-2];
                subdivMesh->HalfEdges[hIndex-2].twin = &subdivMesh->HalfEdges[hIndex+1];
            }

            // For edge points
            subdivMesh->HalfEdges[2*t].target->out = &subdivMesh->HalfEdges[hIndex];

            hIndex += 2;
            fIndex++;
            currentEdge = currentEdge->next;
        }

        subdivMesh->HalfEdges[hIndex-2*n+1].twin = &subdivMesh->HalfEdges[hIndex-2];
        subdivMesh->HalfEdges[hIndex-2].twin = &subdivMesh->HalfEdges[hIndex-2*n+1];

        // For face points
        subdivMesh->Vertices[k].out = &subdivMesh->HalfEdges[hIndex-1];

    }

    qDebug() << " * Created faces and remaining halfedges";

    // For vertex points
    for (k=0; k<numVerts; k++) {
        subdivMesh->Vertices[numFaces + k].out = &subdivMesh->HalfEdges[ 2*inputMesh->Vertices[k].out->index ];
    }

    qDebug() << " * Completed!";
    qDebug() << "   # Vertices:" << subdivMesh->Vertices.size();
    qDebug() << "   # HalfEdges:" << subdivMesh->HalfEdges.size();
    qDebug() << "   # Faces:" << subdivMesh->Faces.size();
//    for(int i = 0; i < subdivMesh->Vertices.size(); i ++){
//        if(subdivMesh->Vertices[i].boundary)
//            qDebug() << "wrong!!!!!";
//    }


}

// Calls setLinearBlendWeights with regular val. 4 and boundary regular val. 3
void CatClarkSubdivider::setLinearBlendWeights(Mesh* mesh, int sharpEV, bool sharpEnd) {
    Subdivider::setLinearBlendWeights(mesh, 4, 3, sharpEV, sharpEnd);
}

// Function that sets the initial blend weights.
void CatClarkSubdivider::setInitialBlendWeights(Mesh* mesh, bool initZeros, short subdivApproach, int sharpEV, bool sharpEnd) {
  //  qDebug() << "set initial blend weights, initZeros =" << initZeros
             //<< "subdivApproach" << subdivApproach << "sharpEV" << sharpEV << "sharpEnd" << sharpEnd ;
    // enum SharpOption {REGULAR, BOUNDARY, ZERO};

    Vertex *v;
    HalfEdge *currentEdge, *faceEdge;
    float fac2, vw, faceval, valfloat, ownLimitWeight;
    bool sharp;

    if (subdivApproach == 3) {
        for(int k = 0; k < mesh->Vertices.size(); k++) {
            v = mesh->Vertices[k].out->twin->target;
            v->origBlendWeight = 1;
        }
        return;
    }

    if (initZeros) {
        for(int k = 0; k < mesh->Vertices.size(); k++) {
            v = mesh->Vertices[k].out->twin->target;
            v->origBlendWeight = 0.0;
        }
    }

    for(int k = 0; k < mesh->Vertices.size(); k++) {
        v = mesh->Vertices[k].out->twin->target;
        valfloat = v->val;
        //sharp or at dart, treat dart as sharp
        sharp = (v->sharp || (v->sharpEnd && sharpEnd));

        //sharp is 0
        if (sharp && sharpEV == 2) { //qDebug() << "sharp is 0=====";
            v->origBlendWeight = 0;
        }
        //sharp is boundary, not correct
        else if (v->boundary ||
                 (sharp && sharpEV == 1)) { // if boundary vertex
            // if irregular vertex: blend weight = 1/(4/6) = 6/4
            v->origBlendWeight = (1.0 - (v->val == 3)) * (6.0/4.0);
            //qDebug() << "sharp is boundary=====?" << v->boundary << "sharp?" << sharp << sharpEV;
        }
        else {//qDebug() << "===========";
            // if irregular vertex
            if(v->val != 4) {

              switch (subdivApproach) {

              //No blend
              case 0:
                // Initialise EV with weight 1 (subdivision approach mentioned/hinted at in original work)
                  //qDebug() << "No blend===========";
                v->origBlendWeight = 1.0;
                break;

                //Linear interpolation
              case 1:
                // Initialise EV with inverse limit stencil. Issue when EVs are not sufficiently separated (multiple EVs per face)
qDebug() << "Linear blend===========";
                // compute limit stencil
                vw = (valfloat-3.0) / (valfloat+5.0);
                fac2 = 4.0 / (valfloat*(valfloat+5.0));
                currentEdge = v->out;
                for(short i = 0; i < v->val; i++) {
                  faceval = currentEdge->polygon->val; // valency of current face
                  vw += fac2 / faceval;
                  currentEdge = currentEdge->twin->next;
                }

                ownLimitWeight = vw + valfloat*fac2/2.0;

                //blend weight = 1/ownLimitWeight
                v->origBlendWeight = 1.0/ownLimitWeight;

                break;

                //Subdiv blend
              case 2:
                // Initialise EV and one-ring NBH to 1. No issues with neighbouring EVs. EFs not directly supported by this!
//qDebug() << "Subdiv blend===========";
                v->origBlendWeight = 1;

                currentEdge = v->out;

                // Set everything in the 1-ring nbh to 1 (not just the star, which would be the loop above...)
                for (int k=0; k<valfloat; k++) {
                  // Vert in the star
                  currentEdge->target->origBlendWeight = 1;

                  faceEdge = currentEdge->next;

                  // Only the verts that are not in the star have to be updated.
                  // So not the vert itself, the star vert, or (after one round) the other star vert. Therefore, -3.
                  for (int l=0; l<faceEdge->polygon->val - 3; l++) {
                    faceEdge->target->origBlendWeight = 1;
                    faceEdge = faceEdge->next;
                  }

                  currentEdge = currentEdge->twin->next;
                }

                break;
              }//switch
            }//if(v->val != 4)
        }
    }
}

// Set control vertices for limit stencil
// [Loop et al., 2009]
// Approximating Subdivision Sur-faces with Gregory Patches for Hardware Tessellation
void CatClarkSubdivider::setLimitSupport(Mesh* mesh) {
    Vertex *v, *v2;
    HalfEdge *currentEdge, *faceEdge;
    Face *face;
    float valfloat, fac1, fac2, vw, faceval;

    for(int k = 0; k < mesh->Vertices.size(); k++) {
        v = mesh->Vertices[k].out->twin->target;
        valfloat = v->val;

        v->limitVertices.clear();
        v->limitWeights.clear();

        // if boundary vertex
        if (HalfEdge* boundaryEdge = vertOnBoundary(v)) {
            // apply general boundary stencil (same for Loop)
            setBoundaryLimitSupport(boundaryEdge, v);
        } else { // inner vertex: limit stencil inner vertex
            v->limitVertices.reserve(1+3*v->val);
            v->limitWeights.reserve(1+3*v->val);

            fac1 = (valfloat-3.0) / (valfloat+5.0);
            vw = fac1;
            fac2 = 4.0 / (valfloat*(valfloat+5.0));
            currentEdge = v->out;
            for(short i = 0; i < v->val; i++) {
                v2 = currentEdge->target;
                face = currentEdge->polygon;
                faceval = face->val;

                // add adjacent vertex
                v->limitVertices.append(v2);
                v->limitWeights.append(fac2/2.0 + fac2/faceval);

                vw += fac2/faceval;

                // Add all incident vertices of face for face point
                faceEdge = currentEdge->next;
                for(short j = 1; j < face->val-1; j++) {
                    v->limitVertices.append(faceEdge->target);
                    v->limitWeights.append(fac2/faceval);
                    faceEdge = faceEdge->next;
                }
                assert(faceEdge->target == v);
                currentEdge = currentEdge->twin->next;
            }

            // add vertex point itself (parent, previous mesh)
            v->limitVertices.append(v);
            v->limitWeights.append(vw + valfloat*fac2/2.0);
        }

#ifdef QT_DEBUG
        float totalweight = 0;
        for(short i = 0; i < v->limitVertices.size(); i++) {
            totalweight += v->limitWeights[i];
        }
        assert(totalweight < 1.01 && totalweight > 0.99);
#endif
    }

    mesh->limitSupportSet = true;
}

// Returns the set of edges connected to the specified vertex that have a positive sharpness value.
QVector<HalfEdge*> CatClarkSubdivider::getSharpEdges(Vertex* vertex) {
    QVector<HalfEdge*> sharpEdges = QVector<HalfEdge*>();
    HalfEdge* currentEdge = vertex->out;

    do {
        if (currentEdge->sharpness > 0) {
            sharpEdges.append(currentEdge);
        }
        currentEdge = currentEdge->twin->next;
    } while (currentEdge != vertex->out);

    return sharpEdges;
}

// Create the set of subdivision masks for the various regions around a sharp vertex point.
void CatClarkSubdivider::setSharpMaskVertexPoint(Vertex* vertex, HalfEdge* firstEdge, Mesh* mesh) {
    QVector<QVector<Vertex*>> sharpMaskVertices(1);
    QVector<QVector<double>> sharpMaskWeights(1);
    QVector<QVector<Face*>> sharpMaskSides(1);
    QVector<QVector<int>> sharpMaskFaces(1);
    HalfEdge* currentEdge = firstEdge;
    Vertex* parent = firstEdge->twin->target;

    sharpMaskVertices.last().append(parent);
    sharpMaskWeights.last().append(6);
    sharpMaskSides.last().append(currentEdge->polygon);

    do {
        if (currentEdge->sharpness > 0) {
            sharpMaskVertices.last().append(currentEdge->target); //Add target vertex to current mask
            sharpMaskWeights.last().append(1);
            sharpMaskSides.last().append(currentEdge->twin->polygon);

            sharpMaskVertices.append(QVector<Vertex*>()); //Create new mask
            sharpMaskWeights.append(QVector<double>());
            sharpMaskFaces.append(QVector<int>());
            sharpMaskSides.append(QVector<Face*>());

            sharpMaskVertices.last().append(parent);
            sharpMaskWeights.last().append(6);
            sharpMaskSides.last().append(currentEdge->polygon);

            sharpMaskVertices.last().append(currentEdge->target); //Add target vertex to new current mask
            sharpMaskWeights.last().append(1);
            sharpMaskSides.last().append(currentEdge->polygon);
        }

        sharpMaskFaces.last().append(currentEdge->polygon->index); //add face to current mask
        currentEdge = currentEdge->prev->twin; //Go to next edge
    } while (currentEdge != firstEdge);

    // Combine first and last mask as they are actually 2 halves of the same mask
    sharpMaskVertices.first().removeFirst(); //remove duplicate parent
    sharpMaskWeights.first().removeFirst();
    sharpMaskSides.first().removeFirst();
    sharpMaskVertices[0].append(sharpMaskVertices.last());
    sharpMaskWeights[0].append(sharpMaskWeights.last());
    sharpMaskFaces[0].append(sharpMaskFaces.last());
    sharpMaskSides[0].append(sharpMaskSides.last());
    sharpMaskVertices.removeLast();
    sharpMaskWeights.removeLast();
    sharpMaskFaces.removeLast();
    sharpMaskSides.removeLast();

    // With more than 2 sharp edges the new vertex point depends solely on its parent
    if (sharpMaskVertices.size() > 2) {
        for (int i = 0; i < sharpMaskVertices.size(); i++) {
            sharpMaskVertices[i] = QVector<Vertex*> {parent};
            sharpMaskWeights[i] = {1};
        }
    }

    // Normalize weights
    for (int i = 0; i < sharpMaskWeights.size(); i++) {
        double sum = 0;
        for (int j = 0; j < sharpMaskWeights[i].size(); j++) {
            sum += sharpMaskWeights[i][j];
        }
        for (int j = 0; j < sharpMaskWeights[i].size(); j++) {
            sharpMaskWeights[i][j] /= sum;
        }
    }

    vertex->setSharpMasks(sharpMaskVertices, sharpMaskWeights, sharpMaskFaces, sharpMaskSides);
}

// Create the set of subdivision masks for the various regions around a sharp edge point.
void CatClarkSubdivider::setSharpMaskEdgePoint(Vertex* vertex, HalfEdge* firstEdge, Mesh* mesh) {
    QVector<QVector<Vertex*>> sharpMaskVertices(2);
    QVector<QVector<double>> sharpMaskWeights(2);
    QVector<QVector<Face*>> sharpMaskSides(2);
    QVector<QVector<int>> sharpMaskFaces(2);
    HalfEdge* currentEdge = firstEdge;

    sharpMaskVertices[0] = QVector<Vertex*> {firstEdge->target, firstEdge->twin->target};
    sharpMaskVertices[1] = QVector<Vertex*> {firstEdge->target, firstEdge->twin->target};
    sharpMaskWeights[0] = QVector<double> {0.5, 0.5};
    sharpMaskWeights[1] = QVector<double> {0.5, 0.5};
    sharpMaskSides[0] = QVector<Face*> {firstEdge->polygon, firstEdge->polygon};
    sharpMaskSides[1] = QVector<Face*> {firstEdge->twin->polygon, firstEdge->twin->polygon};
    sharpMaskFaces[0] = QVector<int> {firstEdge->polygon->index, firstEdge->polygon->index};
    sharpMaskFaces[1] = QVector<int> {firstEdge->twin->polygon->index, firstEdge->twin->polygon->index};

    vertex->setSharpMasks(sharpMaskVertices, sharpMaskWeights, sharpMaskFaces, sharpMaskSides);

//    for (int k = 0; k < sharpMaskVertices.size(); k++) {
////        inputNorms.clear();
//        for(int j = 0; j < sharpMaskVertices[k].size(); j++) {
//            int idx = sharpMaskSides[k][j]->parent;
//            if (idx == -1) {
////                idx = v->sharpMaskSides[k][j]->index;
////                qDebug() << "-1!!!!!!!!!!!!!!!";
//            }// if

////            qDebug() << "idx" << idx <<"sharpMaskSides[k][j]" << sharpMaskSides[k][j]->index;

////            inputNorms.append(v->sharpMaskVertices[k][j]->getSharpNormal(idx));
//        }// for j
////        weights = v->sharpMaskWeights[k];
////        v->sharpLinNormals[k] = linearInterpolation(inputNorms, weights).normalized();
////        wrongNormCnt += (v->linNormal.length() < 0.95);
////        v->sharpSphereNormals[k] = calcSphereNormal(v->sharpLinNormals[k], inputNorms, weights, slerpIterations);
////        v->sharpDotprodLerpSlerps[k] = QVector3D::dotProduct(v->sharpLinNormals[k], v->sharpSphereNormals[k]);
//    }
}

// Function that returns a vertex point given the outgoing edge of parent vertex point
Vertex CatClarkSubdivider::vertexPoint(HalfEdge* firstEdge, Mesh* subdivMesh, short vval, int vindex) {
    short k, n, l;
    HalfEdge* currentEdge;
    Vertex* parentVertex;
    QVector<Vertex *> maskVertices;
    QVector<double> maskWeights;

    QVector<Face*> maskFaces;
    bool sharp = false, willSharp;

    parentVertex = firstEdge->twin->target;
    currentEdge = firstEdge;
    QVector<HalfEdge*> sharpEdges = getSharpEdges(parentVertex);

    float maxSharpness = 0;
    for (int i = 0; i < sharpEdges.size(); i++) {
        maxSharpness = std::max(maxSharpness, sharpEdges[i]->sharpness);
    }
    willSharp = maxSharpness > 1;//qDebug() << "SS" <<willSharp;

    if (sharpEdges.size() == 2) {
        maskVertices.append(sharpEdges[0]->target);
        maskVertices.append(sharpEdges[0]->twin->target);
        maskVertices.append(sharpEdges[1]->target);

        maskFaces.append(sharpEdges[0]->polygon);
        maskFaces.append(sharpEdges[0]->twin->polygon);
        maskFaces.append(sharpEdges[1]->polygon);
        sharp = true;
    }
    else if (HalfEdge* boundaryEdge = vertOnBoundary(parentVertex)) {//qDebug() << "$$$$";
        maskVertices.append(boundaryEdge->target);
        maskVertices.append(boundaryEdge->twin->target);
        maskVertices.append(boundaryEdge->prev->twin->target);
        maskFaces.append(boundaryEdge->polygon);
        maskFaces.append(boundaryEdge->twin->polygon);
        maskFaces.append(boundaryEdge->prev->twin->polygon);
    }

    if (maskVertices.size() > 0) {
        maskWeights.append(1.0/8.0);
        maskWeights.append(6.0/8.0);
        maskWeights.append(1.0/8.0);
        Vertex v(vval, vindex, parentVertex->boundary, maskVertices, maskWeights, maskFaces, 0, parentVertex, sharp);
        if (sharp) {
            setSharpMaskVertexPoint(&v, sharpEdges[0], subdivMesh);
        }
        v.willSharp = willSharp;
        return v;
    }

    QVector<Vertex*> maskverts;
    QVector<Face*> maskfacs;
    int masksize;

    if (sharpEdges.size() > 2) {
        n = sharpEdges.size();
        maskVertices.append(parentVertex);
        maskWeights.append(1);
        maskFaces.append(firstEdge->twin->polygon);
        sharp = true;
        Vertex v(vval, vindex, parentVertex->boundary, maskVertices, maskWeights, maskFaces, 0, parentVertex, sharp);
        setSharpMaskVertexPoint(&v, sharpEdges[0], subdivMesh);
        v.willSharp = willSharp;
        return v;
    }

    n = parentVertex->val;

    for (k=0; k<n; k++) {

        maskVertices.append(currentEdge->target);
        maskWeights.append(1.0/(n*n));
        maskFaces.append(currentEdge->polygon);

        Vertex *faceVertex = &(subdivMesh->Vertices[currentEdge->polygon->index]);
        maskverts = faceVertex->getMaskVertices();
        maskfacs = faceVertex->getMaskFaces();
        masksize = maskverts.size();

        // add support of face-vertex
        for(l = 0; l < masksize; l++) {
            maskVertices.append(maskverts[l]);
            maskFaces.append(maskfacs[l]);
            maskWeights.append(1.0/(n*n*masksize));
        }

        currentEdge = currentEdge->prev->twin;
    }

    maskVertices.append(parentVertex);
    maskWeights.append((n-2.0)/n);
    maskFaces.append(firstEdge->twin->polygon);//qDebug() << "not sharpo" << willSharp;

    Vertex v(vval, vindex, parentVertex->boundary, maskVertices, maskWeights, maskFaces, 0, parentVertex, sharp);
    if (sharp) {
        setSharpMaskVertexPoint(&v, sharpEdges[0], subdivMesh);
    }
    v.willSharp = willSharp;
    if (sharpEdges.size() == 1) {
        v.sharpEnd = true;
    }
//    if(v.boundary)
//        qDebug() << "^^^^^" << v.boundary;
    return v;
}

// Function that returns an edge point, given the edge
Vertex CatClarkSubdivider::edgePoint(HalfEdge* firstEdge, Mesh* subdivMesh, short vval, int vindex) {
    HalfEdge* currentEdge;
    QVector<Vertex *> maskVertices;
    QVector<double> maskWeights;
    QVector<Face*> maskFaces;

    currentEdge = firstEdge;

    // Check if boundary point -> apply general boundary stencil
//    qDebug() << "curent edge on b?" << currentEdge->polygon;
    if (!currentEdge->polygon || !currentEdge->twin->polygon || currentEdge->sharpness > 0) {
        Vertex BEP = boundaryEdgePoint(currentEdge, vval,  vindex);
        BEP.sharp = (currentEdge->sharpness > 0);
        if (BEP.sharp) {
            setSharpMaskEdgePoint(&BEP, currentEdge, subdivMesh);
            BEP.boundary = false;
            BEP.willSharp = currentEdge->sharpness > 1;
//            qDebug() << "7777???BEPPPP" << BEP.boundary;
            return BEP;
        }
        else{
            return BEP;
        }
    }

    Vertex faceVert1, faceVert2;
    int totalMaskSize;
    QVector<Vertex*> maskverts1, maskverts2;
    QVector<Face*> maskFaces1, maskFaces2;
    int masksize1, masksize2, l;
    float w1, w2;

    faceVert1 = subdivMesh->Vertices[currentEdge->polygon->index];
    faceVert2 = subdivMesh->Vertices[currentEdge->twin->polygon->index];

    // get support vertices of face-vertices
    maskverts1 = faceVert1.getMaskVertices();
    maskFaces1 = faceVert1.getMaskFaces();
    masksize1 = maskverts1.size();
    w1 = 1.0 / (4.0*masksize1);
    maskverts2 = faceVert2.getMaskVertices();
    maskFaces2 = faceVert2.getMaskFaces();
    masksize2 = maskverts2.size();
    w2 = 1.0 / (4.0*masksize2);

    totalMaskSize = 2 + masksize1 + masksize2;

    // reserve memory for support vertices
    maskVertices.reserve(totalMaskSize);
    maskWeights.reserve(totalMaskSize);

    // add two vertices from boundary
    maskVertices.append(currentEdge->target);
    maskVertices.append(currentEdge->twin->target);
    maskFaces.append(currentEdge->polygon);
    maskFaces.append(currentEdge->twin->polygon);
    maskWeights.append(1.0/4.0);
    maskWeights.append(1.0/4.0);

    // add support vertices of face-vertex 1
    for(l = 0; l < masksize1; l++) {
        maskVertices.append(maskverts1[l]);
        maskFaces.append(maskFaces1[l]);
        maskWeights.append(w1);
    }

    // add support vertices of face-vertex 2
    for(l = 0; l < masksize2; l++) {
        maskVertices.append(maskverts2[l]);
        maskFaces.append(maskFaces2[l]);
        maskWeights.append(w2);
    }

    Vertex vert = Vertex(vval, vindex, false, maskVertices, maskWeights, maskFaces, 1);
//    if(vert.boundary) {
//        qDebug() << "???BBBBSSS" << vert.boundary;
//        }

    // Set edgepoints as vertices for interpolation (LINEAR) blend function
    vert.linblendParents.reserve(2);
    vert.linblendParents.append(firstEdge->target);
    vert.linblendParents.append(firstEdge->twin->target);
//qDebug() << "7777???BBBBSSS" << vert.boundary;
    return vert;
}

// Function that returns an edge point, given one of the inner edges
Vertex CatClarkSubdivider::facePoint(HalfEdge* firstEdge, short vval, int vindex) {
    short k, n;
    HalfEdge* currentEdge;
    QVector<Vertex *> maskVertices;
    QVector<double> maskWeights;
    QVector<Face*> maskFaces;
    double weight;
    Vertex *v;

    n = firstEdge->polygon->val;
    weight = 1.0/n;

    currentEdge = firstEdge;

    maskVertices.reserve(n);
    maskWeights.reserve(n);

    // Add all incident vertices of face as control vertices
    // by visiting the inner edges of the face
    for (k = 0; k < n; k++) {
        v = currentEdge->target;
        maskVertices.append(v);
        maskWeights.append(weight);
        maskFaces.append(currentEdge->polygon);

        currentEdge = currentEdge->next;
    }

    Vertex vert = Vertex(vval, vindex, false, maskVertices, maskWeights, maskFaces, 2);

    // Set incident vertices for interpolation (LINEAR) blend function
    vert.linblendParents = maskVertices;
    if(vert.boundary)
        qDebug() << "******" << vert.boundary;
    return vert;
}
