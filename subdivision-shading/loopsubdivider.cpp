#include "loopsubdivider.h"
#include "math.h"
#include "assert.h"


LoopSubdivider::LoopSubdivider() {
    qDebug() << "✓✓ LoopSubdivider constructor (Empty)";
}

// Implementation of Subdivider::subdivide() using the Loop
// subdivision scheme
void LoopSubdivider::subdivide(Mesh* inputMesh, Mesh* subdivMesh) {
    subdivMesh->level = inputMesh->level + 1;

    int numVerts, numHalfEdges, numFaces;
    int k;

    qDebug() << ":: Creating new Loop mesh";

    numVerts = inputMesh->Vertices.size();
    numHalfEdges = inputMesh->HalfEdges.size();
    numFaces = inputMesh->Faces.size();

    // Reserve memory
    subdivMesh->Vertices.reserve(numVerts + numHalfEdges/2);
    subdivMesh->HalfEdges.reserve(2*numHalfEdges + 6*numFaces);
    subdivMesh->Faces.reserve(4*numFaces);

    // Create vertex points
    createVertexPoints(inputMesh, subdivMesh);

    qDebug() << " * Created vertex points";

    // Create edge points
    qDebug() << " * Creating edge points";
    createEdgePoints(inputMesh, subdivMesh, numHalfEdges);

    qDebug() << " * Created edge points";

    // Split halfedges
    splitHalfEdges(inputMesh, subdivMesh, numHalfEdges, numVerts, 0);

    qDebug() << " * Split halfedges";


    // Create faces and remaining halfedges
    createFaces(inputMesh, subdivMesh, numHalfEdges, numFaces);

    qDebug() << " * Created faces";

    // For vertex points
    for (k=0; k<numVerts; k++) {
        subdivMesh->Vertices[k].out = &subdivMesh->HalfEdges[ 2*inputMesh->Vertices[k].out->index ];
    }

    qDebug() << " * Completed!";
    qDebug() << "   # Vertices:" << subdivMesh->Vertices.size();
    qDebug() << "   # HalfEdges:" << subdivMesh->HalfEdges.size();
    qDebug() << "   # Faces:" << subdivMesh->Faces.size();
}

// Function that returns a vertex point given the outgoing edge of parent vertex point
Vertex LoopSubdivider::vertexPoint(HalfEdge* firstEdge, short vval, int vindex) {
    float stencilValue, stencilPart;
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
    willSharp = maxSharpness > 1;

    if (sharpEdges.size() == 2) {
        maskVertices.append(sharpEdges[0]->target);
        maskVertices.append(sharpEdges[0]->twin->target);
        maskVertices.append(sharpEdges[1]->target);

        maskFaces.append(sharpEdges[0]->polygon);
        maskFaces.append(sharpEdges[0]->twin->polygon);
        maskFaces.append(sharpEdges[1]->polygon);
        sharp = true;
    } else if (HalfEdge* boundaryEdge = vertOnBoundary(parentVertex)) {
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
            setSharpMaskVertexPoint(&v, sharpEdges[0]);
        }
        v.willSharp = willSharp;
        return v;
    }

    if (sharpEdges.size() > 2) {
        int n = sharpEdges.size();
        maskVertices.append(parentVertex);
        maskWeights.append(1);
        maskFaces.append(firstEdge->twin->polygon);
        sharp = true;
        Vertex v(vval, vindex, parentVertex->boundary, maskVertices, maskWeights, maskFaces, 0, parentVertex, sharp);
        setSharpMaskVertexPoint(&v, sharpEdges[0]);
        v.willSharp = willSharp;
        return v;
    }

    // Check if boundary point -> apply general boundary stencil
    if (HalfEdge* boundaryEdge = vertOnBoundary(parentVertex)) {
//        return boundaryVertexPoint(boundaryEdge, parentVertex, vval, vindex);
    }

    // LOOPS RULES
    stencilPart = 3.0/8.0 + (1.0/4.0)*cos(2.0*M_PI / vval);
    stencilValue = (1.0 / vval) * (5.0/8.0 - stencilPart*stencilPart);

    maskVertices.reserve(vval);
    maskWeights.reserve(vval);
    maskFaces.reserve(vval);

    // parent vertex point (previous mesh)
    maskVertices.append(parentVertex);
    maskWeights.append(1.0 - vval*stencilValue);
    maskFaces.append(firstEdge -> twin -> polygon);

//    int n = parentVertex->val;

    for (short k=0; k<vval; k++) {
        // add support of adjacent vertices
        maskVertices.append(currentEdge->target);
        maskWeights.append(stencilValue);
        maskFaces.append(currentEdge -> polygon);

        currentEdge = currentEdge->twin->next;
    }

    Vertex v(vval, vindex, parentVertex->boundary, maskVertices, maskWeights, maskFaces, 0, parentVertex);
    if (sharp) {
        setSharpMaskVertexPoint(&v, sharpEdges[0]);
    }
    v.willSharp = willSharp;
    if (sharpEdges.size() == 1) {
        v.sharpEnd = true;
    }
    return v;

}

// Function that returns an edge point, given the edge
Vertex LoopSubdivider::edgePoint(HalfEdge* firstEdge, short vval, int vindex) {
    QVector<Vertex *> maskVertices;
    QVector<double> maskWeights;
    QVector<Face*> maskFaces;
    HalfEdge* currentEdge;

    currentEdge = firstEdge;

    // Check if boundary point -> apply general boundary stencil
    if(vval == 4 || currentEdge->sharpness > 0) { //qDebug() << "vval == 4 or sharp?";
        Vertex BEP = boundaryEdgePoint(currentEdge, vval, vindex);
//        return boundaryEdgePoint(firstEdge, vval, vindex);
        BEP.sharp = currentEdge->sharpness > 0;
        if(BEP.sharp) {
            setSharpMaskEdgePoint(&BEP, currentEdge);
            BEP.boundary = false;
            BEP.willSharp = currentEdge -> sharpness > 1;
            return BEP;
        }
        else{ //qDebug() << "boundary##";
            return BEP;
        }

    }

    // Apply Loop edge stencil
    maskVertices.reserve(4);
    maskWeights.reserve(4);
    maskFaces.reserve(4);

    maskVertices.append(firstEdge->target); // edge endpoint 1
    maskVertices.append(firstEdge->next->target);
    maskVertices.append(firstEdge->twin->target); // edge endpoint 2
    maskVertices.append(firstEdge->twin->next->target);

    maskWeights.append(6.0/16.0);
    maskWeights.append(2.0/16.0);
    maskWeights.append(6.0/16.0);
    maskWeights.append(2.0/16.0);

    maskFaces.append(currentEdge -> polygon);
    maskFaces.append(currentEdge -> next -> polygon);
    maskFaces.append(currentEdge -> twin -> polygon);
    maskFaces.append(currentEdge -> twin -> next -> polygon);

    Vertex vert = Vertex(vval, vindex, false, maskVertices, maskWeights, maskFaces, 1);

    // add support for LINEAR BLEND function
    vert.linblendParents.reserve(2);
    vert.linblendParents.append(firstEdge->target);
    vert.linblendParents.append(firstEdge->twin->target);

    return vert;
}

// Function that adds all vertex points to subdivMesh, given previous (inputMesh)
void LoopSubdivider::createVertexPoints(Mesh* inputMesh, Mesh* subdivMesh) {
    short n;

    qDebug() << ":: Creating vertex points";

    // Create vertex points
    for (int k=0; k<inputMesh->Vertices.size(); k++) {
        Vertex oldVertex = inputMesh->Vertices[k];
        n = oldVertex.val;

        // Out, Valence, Index
        Vertex newVertexPoint = vertexPoint(oldVertex.out, n, k);
        subdivMesh->Vertices.append(newVertexPoint);

        // Add vertex point to its predecessors
        // Needed when normal subdivision > level then geometry subdivision
        subdivMesh->Vertices[k].addToParents();
    }
}

// Function that adds all edge points to subdivMesh, given previous (inputMesh)
void LoopSubdivider::createEdgePoints(Mesh* inputMesh, Mesh* subdivMesh, int numHalfEdges) {
    int k, val;
    int vIndex;
    HalfEdge* currentEdge;

    vIndex = inputMesh->Vertices.size();

    // Create edge points
    for (k=0; k < numHalfEdges; k++) {
        currentEdge = &inputMesh->HalfEdges[k];

        if (k < currentEdge->twin->index) {
            // Out, Valence, Index
            // If on boundary -> valency = 4
            val = ((!currentEdge->polygon || !currentEdge->twin->polygon) ? 4 : 6);

            subdivMesh->Vertices.append( edgePoint(currentEdge, val, vIndex) );

            vIndex++;
        }//if
    }//for
}

// Function that creates all faces for subdivMesh, given previous (inputMesh)
void LoopSubdivider::createFaces(Mesh* inputMesh, Mesh* subdivMesh, int numHalfEdges, int numFaces) {
    int k, m, s, t;
    int hIndex, fIndex;
    HalfEdge* currentEdge;

    hIndex = 2*numHalfEdges;
    fIndex = 0;

    for (k=0; k<numFaces; k++) {
        currentEdge = inputMesh->Faces[k].side;

        // Three outer faces

        for (m=0; m<3; m++) {

            s = currentEdge->prev->index;
            t = currentEdge->index;

            // Side, Val, Index
            subdivMesh->Faces.append(Face(3, fIndex, inputMesh->Faces[k].index));

            subdivMesh->Faces[fIndex].side = &subdivMesh->HalfEdges[ 2*t ];

            // Target, Next, Prev, Twin, Poly, Index
            subdivMesh->HalfEdges.append(HalfEdge( subdivMesh->HalfEdges[2*s].target,
                                         &subdivMesh->HalfEdges[2*s+1],
                    &subdivMesh->HalfEdges[ 2*t ],
                    nullptr,
                    &subdivMesh->Faces[fIndex],
                    hIndex ));

            subdivMesh->HalfEdges.append(HalfEdge( subdivMesh->HalfEdges[2*t].target,
                                         nullptr,
                                         nullptr,
                                         &subdivMesh->HalfEdges[hIndex],
                                         nullptr,
                                         hIndex+1 ));

            subdivMesh->HalfEdges[hIndex].twin = &subdivMesh->HalfEdges[hIndex+1];

            subdivMesh->HalfEdges[2*s+1].next = &subdivMesh->HalfEdges[2*t];
            subdivMesh->HalfEdges[2*s+1].prev = &subdivMesh->HalfEdges[hIndex];
            subdivMesh->HalfEdges[2*s+1].polygon = &subdivMesh->Faces[fIndex];

            subdivMesh->HalfEdges[2*t].next = &subdivMesh->HalfEdges[hIndex];
            subdivMesh->HalfEdges[2*t].prev = &subdivMesh->HalfEdges[2*s+1];
            subdivMesh->HalfEdges[2*t].polygon = &subdivMesh->Faces[fIndex];

            // For edge points
            subdivMesh->HalfEdges[2*t].target->out = &subdivMesh->HalfEdges[hIndex];

            hIndex += 2;
            fIndex++;
            currentEdge = currentEdge->next;

        }

        // Inner face

        // Side, Val, Index
        subdivMesh->Faces.append(Face(3, fIndex, inputMesh->Faces[k].index));

        subdivMesh->Faces[fIndex].side = &subdivMesh->HalfEdges[ hIndex-1 ];

        for (m=0; m<3; m++) {
            if (m==2) {
                subdivMesh->HalfEdges[hIndex - 1].next = &subdivMesh->HalfEdges[hIndex - 5];
            } else {
                subdivMesh->HalfEdges[hIndex - 5 + 2*m].next = &subdivMesh->HalfEdges[hIndex - 5 + 2*(m+1)];
            }

            if (m==0) {
                subdivMesh->HalfEdges[hIndex - 5].prev = &subdivMesh->HalfEdges[hIndex - 1];
            } else {
                subdivMesh->HalfEdges[hIndex - 5 + 2*m].prev = &subdivMesh->HalfEdges[hIndex - 5 + 2*(m-1)];
            }
            subdivMesh->HalfEdges[hIndex - 5 + 2*m].polygon = &subdivMesh->Faces[fIndex];
        }
        fIndex++;
    }
}

void LoopSubdivider::setLinearBlendWeights(Mesh* mesh, int sharpEV, bool sharpEnd) {
    // Call setLinearBlendWeights with regular val. 6 and boundary regular val. 4
    Subdivider::setLinearBlendWeights(mesh, 6, 4, sharpEV, sharpEnd);
}


void LoopSubdivider::setInitialBlendWeights(Mesh* mesh, bool initZeros, short subdivApproach, int sharpEV, bool sharpEnd) {
    //qDebug() << "set initial blend weights, initZeros =" << initZeros
                 //<< "subdivApproach" << subdivApproach << "sharpEV" << sharpEV << "sharpEnd" << sharpEnd ;

    Vertex *v;
    HalfEdge *currentEdge, *faceEdge;
    float n, wn, wnpart2, divisor, valfloat;
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


//        if(v->boundary && v->val != 4) {
//            v->origBlendWeight = 6.0 / 4.0;
//        }

        //sharp is 0
        if (sharp && sharpEV == 2) { //qDebug() << "sharp is 0=====";
            v->origBlendWeight = 0;
        }
        //sharp is boundary, not correct
        else if (v->boundary ||
         (sharp && sharpEV == 1)) { // if boundary vertex
            // if irregular vertex: blend weight = 1/(4/6) = 6/4
            v->origBlendWeight = (1.0 - (v->val == 4)) * (6.0/4.0);
//            v->origBlendWeight = 6.0 / 4.0;
            //qDebug() << "sharp is boundary=====" << v->boundary << sharp << sharpEV;
        }

        //sharp is boundary
        else if (v->boundary ||
                 (sharp && sharpEV == 1)) { // if boundary vertex
            // if irregular vertex: blend weight = 1/(4/6) = 6/4
            v->origBlendWeight = (1.0 - (v->val == 3)) * (6.0/4.0);
            //qDebug() << "sharp is boundary=====?" << v->boundary << "sharp?" << sharp << sharpEV;
        }

        //irregular
        else if(!v->boundary && v->val != 6) {
            n = v->val;

            //            qDebug() << "orig blendweight index " << v->index << " with valency:" << v->val << ":" << v->origBlendWeight;

        switch (subdivApproach) {

        case 0:
          // Initialise EV with weight 1 (subdivision approach mentioned/hinted at in original work)
          v->origBlendWeight = 1.0;
          break;

        case 1:
          // Initialise EV with inverse limit stencil. Issue when EVs are not sufficiently separated (multiple EVs per face)

          // compute limit stencil
            wnpart2 = (3.0/8.0) + (2.0/8.0)*cos(2.0*M_PI/n);
            wn = (3.0*n)/(5.0 - 8.0*wnpart2*wnpart2);
            divisor = n+wn;
            v->origBlendWeight = divisor / wn;
          break;

        case 2:
          // Initialise EV and one-ring NBH to 1. No issues with neighbouring EVs. EFs not directly supported by this!

          v->origBlendWeight = 1;

          currentEdge = v->out;

          // Set everything in the 1-ring nbh to 1 (not just the star, which would be the loop above...)
          for (int k=0; k<n; k++) {
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
        }
      } else {
          //v->origBlendWeight = 0.0;
          //Nothing.
      }

    }

}


// Set control vertices for limit stencil
// [Li et al., 2011] Approximation of Loop subdivision surfacesfor fast rendering
void LoopSubdivider::setLimitSupport(Mesh* mesh) {
    Vertex *v, *v2;
    HalfEdge *currentEdge;
    float valfloat, fac1, fac2, wn, wnpart2, divisor;
    for(int k = 0; k < mesh->Vertices.size(); k++) {
        v = mesh->Vertices[k].out->twin->target;

        v->limitVertices.clear();
        v->limitWeights.clear();

        // if boundary vertex
        if (HalfEdge* boundaryEdge = vertOnBoundary(v)) {
            // apply general boundary stencil (same for Catmull-Clark)
            setBoundaryLimitSupport(boundaryEdge, v);
        } else { // Inner vertex: compute limit stencil for Loop inner vertex
            valfloat = v->val;
            wnpart2 = (3.0/8.0) + ((2.0/8.0)*cos(2.0*M_PI/valfloat));
            wn = (3.0*valfloat)/(5.0 - (8.0*wnpart2*wnpart2));
            divisor = valfloat+wn;
            fac1 = wn/divisor;
            fac2 = 1/divisor;

            v->limitVertices.reserve(1+v->val);
            v->limitWeights.reserve(1+v->val);
            v->limitVertices.append(v);
            v->limitWeights.append(fac1);

            currentEdge = v->out;
            for(short i = 0; i < v->val; i++) {
                v2 = currentEdge->target;

                v->limitVertices.append(v2);
                v->limitWeights.append(fac2);

                currentEdge = currentEdge->twin->next;
            }
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

// Create the set of subdivision masks for the various regions around a sharp edge point.
void LoopSubdivider::setSharpMaskEdgePoint(Vertex* vertex, HalfEdge* firstEdge) {
    QVector<QVector<Vertex*>> sharpMaskVertices(2);
    QVector<QVector<double>> sharpMaskWeights(2);
    QVector<QVector<Face*>> sharpMaskSides(2);
    QVector<QVector<int>> sharpMaskFaces(2);
    HalfEdge* currentEdge = firstEdge;

    sharpMaskVertices[0] = QVector<Vertex*> {currentEdge->target, currentEdge->twin->target};
    sharpMaskVertices[1] = QVector<Vertex*> {currentEdge->target, currentEdge->twin->target};
    //qDebug() << "sharpMaskVertices[0][0]" << sharpMaskVertices[0][0]->index << "sharpMaskVertices[0][1]" << sharpMaskVertices[0][1]->index;
    //qDebug() << "sharpMaskVertices[1][0]" << sharpMaskVertices[1][0]->index << "sharpMaskVertices[1][1]" << sharpMaskVertices[1][1]->index;

    sharpMaskWeights[0] = QVector<double> {0.5, 0.5};
    sharpMaskWeights[1] = QVector<double> {0.5, 0.5};

    sharpMaskSides[0] = QVector<Face*> {currentEdge->polygon, currentEdge->polygon};
    sharpMaskSides[1] = QVector<Face*> {currentEdge->twin->polygon, currentEdge->twin->polygon};

    sharpMaskFaces[0] = QVector<int> {static_cast<int>(currentEdge->polygon->index), static_cast<int>(currentEdge->polygon->index)};
    sharpMaskFaces[1] = QVector<int> {static_cast<int>(currentEdge->twin->polygon->index), static_cast<int>(currentEdge->twin->polygon->index)};
//        sharpMaskFaces[0] = QVector<int> {currentEdge->polygon->index, currentEdge->polygon->index};
//        sharpMaskFaces[1] = QVector<int> {currentEdge->twin->polygon->index, currentEdge->twin->polygon->index};

    //    sharpMaskFaces[0].append(currentEdge->polygon->index);
    //    sharpMaskFaces[1].append(currentEdge->twin->polygon->index);

    vertex->setSharpMasks(sharpMaskVertices, sharpMaskWeights, sharpMaskFaces, sharpMaskSides);



}

// Create the set of subdivision masks for the various regions around a sharp vertex point.
void LoopSubdivider::setSharpMaskVertexPoint(Vertex* vertex, HalfEdge* firstEdge) {
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

//        sharpMaskFaces.last().append(currentEdge->prev->twin->polygon->index); //add face to current mask
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

// Returns the set of edges connected to the specified vertex that have a positive sharpness value.
QVector<HalfEdge*> LoopSubdivider::getSharpEdges(Vertex* vertex) {
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
