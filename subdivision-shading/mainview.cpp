#include "mainview.h"
#include "math.h"
#include "assert.h"


MainView::MainView(QWidget *Parent) : QOpenGLWidget(Parent) {
    qDebug() << "✓✓ MainView constructor";

    // set defaults
    modelLoaded = false;
    drawMode = WIREFRAME;
    colourMode = REFLECTION;
    isophotesFrequency = 1;
    blendNormals = false;

    rotAngleY = 0.0;
    rotAngleX = 0.0;
    scale = 1.1;
    maxDiff = 1.0;
    dispRatio = 1.0;
    subdivType = Mesh::SubdivType::CATCLARK;
}

MainView::~MainView() {
    qDebug() << "✗✗ MainView destructor";

    glDeleteBuffers(1, &meshCoordsBO);
    glDeleteBuffers(1, &meshNormalsBO);
    glDeleteBuffers(1, &meshColoursBO);
    glDeleteBuffers(1, &meshIndexBO);
    glDeleteVertexArrays(1, &meshVAO);

    debugLogger->stopLogging();

    delete debugLogger;

    delete mainShaderProg;

    delete vertReflectionShader;
    delete reflectionFragShader;
    delete isophotesFragShader;
    delete vertshader;
    delete fixedNormalFragShader;
}

// Compile shader code
void setShaderSource(QOpenGLShader *shader, const char *path) {
    if(shader->compileSourceFile(path)) {
        qDebug() << "Compiled shader with and source path" << path << "successfully";
    } else {
        exit(-1);
    }
}

// Create a new vertex shader
QOpenGLShader *createVertexShader(const char *path) {
    QOpenGLShader *shader = new QOpenGLShader(QOpenGLShader::Vertex);
    setShaderSource(shader, path);
    return shader;
}

// Create a new fragment shader
QOpenGLShader *createFragmentShader(const char *path) {
    QOpenGLShader *shader = new QOpenGLShader(QOpenGLShader::Fragment);
    setShaderSource(shader, path);
    return shader;
}


// Compile all available/needed shaders
void MainView::createShaderPrograms() {
    qDebug() << ".. createShaderPrograms";

    mainShaderProg = new QOpenGLShaderProgram();

//    vertReflectionShader = createVertexShader(":/shaders/vertshader-reflection.glsl");
    reflectionFragShader = createFragmentShader(":/shaders/fragshader-reflection.glsl");
    isophotesFragShader = createFragmentShader(":/shaders/fragshader-isophotes.glsl");
    vertshader = createVertexShader(":/shaders/vertshader.glsl");
    fixedNormalFragShader = createFragmentShader(":/shaders/fragshader-fixednormal.glsl");
}

// Add all shaders needed for current colour mode
void MainView::addShaders() {
    mainShaderProg->removeAllShaders();

    mainShaderProg->addShader(vertshader);

    switch(colourMode) {
    case REFLECTION: // FIXED COLOUR + REFLECTION
        mainShaderProg->addShader(reflectionFragShader);
        break;
    case ISOPHOTES:
        mainShaderProg->addShader(isophotesFragShader);
        break;
    default:
        mainShaderProg->addShader(fixedNormalFragShader);
        break;
    }

    mainShaderProg->link();

    // Create shader uniforms
    uniModelViewMatrix = glGetUniformLocation(mainShaderProg->programId(), "modelviewmatrix");
    uniProjectionMatrix = glGetUniformLocation(mainShaderProg->programId(), "projectionmatrix");
    uniNormalMatrix = glGetUniformLocation(mainShaderProg->programId(), "normalmatrix");

    uniIsoFreq = glGetUniformLocation(mainShaderProg->programId(), "uniIsoFreq");

    uniShadingType = glGetUniformLocation(mainShaderProg->programId(), "shadingType");
}

void MainView::createBuffers() {

    qDebug() << ".. createBuffers";

    glGenVertexArrays(1, &meshVAO);
    glBindVertexArray(meshVAO);

    glGenBuffers(1, &meshCoordsBO);
    glBindBuffer(GL_ARRAY_BUFFER, meshCoordsBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glGenBuffers(1, &meshNormalsBO);
    glBindBuffer(GL_ARRAY_BUFFER, meshNormalsBO);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glGenBuffers(1, &meshIndexBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, meshIndexBO);

    glBindVertexArray(0);
}

// Set the right buffers for the current colour mode
void MainView::setBuffers() {
    glBindBuffer(GL_ARRAY_BUFFER, meshCoordsBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(QVector3D)*vertexCoords2.size(), vertexCoords2.data(), GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, meshNormalsBO);

    if(colourMode == REFLECTION || colourMode == ISOPHOTES) { // FIXED COLOUR + REFLECTION
        glBufferData(GL_ARRAY_BUFFER, sizeof(QVector3D)*vertexNormals2.size(), vertexNormals2.data(), GL_DYNAMIC_DRAW);
    } else { // NORMAL BUFFER or DIFF colours
        glBufferData(GL_ARRAY_BUFFER, sizeof(QVector3D)*vertexColours2.size(), vertexColours2.data(), GL_DYNAMIC_DRAW);
    }

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, meshIndexBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*polyIndices2.size(), polyIndices2.data(), GL_DYNAMIC_DRAW);
}

// Compute color by value from a heat map specified by min, max
QVector3D heatmap(float min, float max, float value) {
    float ratio = fmin(1.0, (value-min) / (max - min));
    float r = ratio;
    float g = 1 - ratio;
    return QVector3D(r, 0, g);
}

// Set the polygon indices for the buffers
void MainView::setPolyIndices(Mesh *mesh) {
    HalfEdge* currentEdge;
    int k, n, m;

    polyIndices.clear();
    polyIndices.reserve(mesh->HalfEdges.size() + mesh->Faces.size());

    for (k=0; k<mesh->Faces.size(); k++) {
        n = mesh->Faces[k].val;
        currentEdge = mesh->Faces[k].side;
        for (m=0; m<n; m++) {
            polyIndices.append(currentEdge->target->index);
            currentEdge = currentEdge->next;
        }
        polyIndices.append(maxInt);
    }

    polyIndices2.clear();
    polyIndices2.reserve(mesh->HalfEdges.size() + mesh->Faces.size());
    int cnt = 0;

    for (k=0; k<mesh->Faces.size(); k++) {
        n = mesh->Faces[k].val;
        for (m=0; m<n; m++) {
            polyIndices2.append(cnt);
            cnt++;
        }
        polyIndices2.append(maxInt);
    }

    meshIBOSize = polyIndices.size();
    meshIBOSize2 = polyIndices2.size();
}

// Copy vertex coords from mesh
void MainView::setCoords(Mesh *mesh) {
    QVector<Vertex> verts = mesh->Vertices;
    int k, n, m;

    vertexCoords.clear();
    vertexCoords.reserve(verts.size());

    for (k=0; k<verts.size(); k++) {
        vertexCoords.append(verts[k].coords);
    }


    HalfEdge* currentEdge;

    vertexCoords2.clear();
    vertexCoords2.reserve(mesh->HalfEdges.size());

    for (k=0; k<mesh->Faces.size(); k++) {
        n = mesh->Faces[k].val;
        currentEdge = mesh->Faces[k].side;
        for (m=0; m<n; m++) {
            vertexCoords2.append(verts[currentEdge->target->index].coords);
            currentEdge = currentEdge->next;
        }
    }
}

// Copy vertex normals from mesh
void MainView::setNormals(Mesh *mesh) {
    QVector<Vertex> verts = mesh->Vertices;
    Vertex v;

    vertexNormals.clear();
    vertexNormals.reserve(verts.size());
    for (int k=0; k<verts.size(); k++) {
        vertexNormals.append(verts[k].normal);
        v = verts[k];
    }

    HalfEdge* currentEdge;
    int k, n, m, par;

    vertexNormals2.clear();
    vertexNormals2.reserve(mesh->HalfEdges.size());

    for (k=0; k<mesh->Faces.size(); k++) {
        n = mesh->Faces[k].val;
        currentEdge = mesh->Faces[k].side;
        for (m=0; m<n; m++) {
            v = verts[currentEdge->target->index];
            par = mesh->Faces[k].parent;
            if (par < 0) {
                par = mesh->Faces[k].index;
            }
            vertexNormals2.append(v.getSharpNormal(par));
            currentEdge = currentEdge->next;
        }
    }
}

// Set max difference value for scaling difference colors
void MainView::setMaxDiff(float value) {
    maxDiff = value;
}

// Copy differences between LERP and SLERP normals
void MainView::setLerpSlerpDiffs(Mesh *mesh) {
    QVector<Vertex> verts = mesh->Vertices;

    lerpSlerpDiffs.clear();
    lerpSlerpDiffs.reserve(verts.size());

    for (int k=0; k<verts.size(); k++) {
        lerpSlerpDiffs.append(verts[k].dotprodLerpSlerp);
    }
}

// Compute differences between surface and subdivided normals
void MainView::setSurfaceSubdivDiffs(Mesh *mesh) {
    QVector<Vertex> verts = mesh->Vertices;

    faceSlerpDiffs.clear();
    faceSlerpDiffs.reserve(verts.size());

    for (int k=0; k<verts.size(); k++) {
        faceSlerpDiffs.append(QVector3D::dotProduct(verts[k].surfaceNormal, verts[k].normal));
    }
}

// Compute differences between left and right view normals
void MainView::setLeftRightDiffs(Mesh *mesh) {
    QVector<Vertex> verts = mesh->Vertices;

    leftRightDiffs.clear();
    leftRightDiffs.reserve(verts.size());

    for (int k=0; k<verts.size(); k++) {
        leftRightDiffs.append(QVector3D::dotProduct(verts[k].leftNormal, verts[k].rightNormal));
    }
}

// Copy differences between normal and normal projected to (subdiv) limit
void MainView::setLimitDiffs(Mesh *mesh) {
    QVector<Vertex> verts = mesh->Vertices;
    int m, n;
    HalfEdge* currentEdge;
    Vertex v;

    limitDiffs.clear();
    limitDiffs.reserve(vertCount);
    limitDiffs2.clear();
    limitDiffs2.reserve(vertCount);

    for (int k=0; k<verts.size(); k++) {
        limitDiffs.append(verts[k].dotprodLimit);
    }

    for (int k=0; k<mesh->Faces.size(); k++) {
        n = mesh->Faces[k].val;
        currentEdge = mesh->Faces[k].side;
        for (m=0; m<n; m++) {
            v = verts[currentEdge->target->index];
            limitDiffs2.append(v.dotprodLimit);
            currentEdge = currentEdge->next;
        }
    }
}

// Copy blend weights
void MainView::setBlendWeights(Mesh *mesh) {
    QVector<Vertex> verts = mesh->Vertices;
    int m, n, par;
    HalfEdge* currentEdge;
    Vertex v;

    blendWeights.clear();
    blendWeights.reserve(vertCount);
    blendWeights2.clear();
    blendWeights2.reserve(vertCount);

    for (int k=0; k<verts.size(); k++) {
        blendWeights.append(verts[k].limitBlendWeight);
    }

    for (int k=0; k<mesh->Faces.size(); k++) {
        n = mesh->Faces[k].val;
        currentEdge = mesh->Faces[k].side;
        for (m=0; m<n; m++) {
            v = verts[currentEdge->target->index];
            par = mesh->Faces[k].parent;
            if (par < 0) {
                par = mesh->Faces[k].index;
            }
            //blendWeights2.append(v.getSharpBlendWeight(par));
            blendWeights2.append(v.limitBlendWeight);
            currentEdge = currentEdge->next;
        }
    }
}

// Check if vertices are regular and save results
void MainView::setVertsRegular(Mesh *mesh) {
    QVector<Vertex> verts = mesh->Vertices;

    int normval;
    int extrordval;

    vertIsRegular.clear();
    vertIsRegular.reserve(verts.size());

    if(subdivType == Mesh::SubdivType::CATCLARK) {
        normval = 4;
        extrordval = 3;
    } else {
        normval = 6;
        extrordval = 4;
    }

    for (int k=0; k<verts.size(); k++) {
        Vertex v = verts[k];
        vertIsRegular.append((!v.boundary && v.val == normval) || (v.boundary && v.val == extrordval));
    }
}

// Compute and set colors by differences between LERP and SLERP normals
void MainView::setColoursByDiffLerpSlerp() {
    for (int k=0; k<vertCount; k++) {
        vertexColours.append(heatmap(0.0,maxDiff,1.0-lerpSlerpDiffs[k]));
    }
}

// Compute and set colors by differences between surface and subdivision normals
void MainView::setColoursByDiffSurfaceSubdiv() {
    for (int k=0; k<vertCount; k++) {
        vertexColours.append(heatmap(0.0,maxDiff,1.0-faceSlerpDiffs[k]));
    }
}

// Compute and set colors by differences between normals of left and right view
void MainView::setColoursByDiffLeftRight() {
    for (int k=0; k<vertCount; k++) {
        vertexColours.append(heatmap(0.0,maxDiff,1.0-leftRightDiffs[k]));
    }
}

// Compute and set colors by differences between normals and their limit projection
void MainView::setColoursByDiffLimit() {
    for (int k=0; k<vertCount; k++) {
        vertexColours.append(heatmap(0.0,maxDiff,1.0-limitDiffs[k]));
    }
    for (int k = 0; k < limitDiffs2.size(); k++) {
        vertexColours2.append(heatmap(0.0,maxDiff,1.0-limitDiffs2[k]));
    }
}

// Compute and set colors by blend weights
void MainView::setColoursByBlendWeights() {
    for (int k=0; k<vertCount; k++) {
        vertexColours.append(heatmap(0.0,1.0,blendWeights[k]));
    }
    for (int k=0; k<blendWeights2.size(); k++) {
        vertexColours2.append(heatmap(0.0,1.0,blendWeights2[k]));
    }
}

// Compute and set colors by normals (normal buffer)
void MainView::setColoursByNormalBuffer() {
    const QVector3D plus = QVector3D(1.2,1.0,1.2);
    for (int k=0; k<vertCount; k++) {
        vertexColours.append((vertexNormals[k] + plus)/2.0);
    }
    for (int k=0; k<vertexNormals2.size(); k++) {
        vertexColours2.append((vertexNormals2[k] + plus)/2.0);
    }
}

// Compute and set colors based on regularity.
// If vertex is regular, color is green
// If vertex is irregular, color is red
void MainView::setColoursByExtrordVerts() {
    for (int k=0; k<vertCount; k++) {

        if(vertIsRegular[k]) {
            vertexColours.append(QVector3D(0.0, 1.0, 0.0));
        } else {
            vertexColours.append(QVector3D(1.0, 0.0, 0.0));
        }
    }
}

// Set colors by given colour mode
void MainView::setColours() {
    vertexColours.clear();
    vertexColours.reserve(vertCount);
    vertexColours2.clear();
    vertexColours2.reserve(vertexNormals2.size());

    // NORMAL COLORS
    switch(colourMode) {
    case NORMALSDIFF:
        setColoursByDiffLerpSlerp();
        break;
    case NORMALSDIFF_MESH:
        setColoursByDiffSurfaceSubdiv();
        break;
    case LIMIT_NORMALS_DIFF:
        setColoursByDiffLimit();
        break;
    case EXTRORDVERTS:
        setColoursByExtrordVerts();
        break;
    case LEFT_RIGHT_DIFF:
        setColoursByDiffLeftRight();
        break;
    case BLENDWEIGHTS:
        setColoursByBlendWeights();
        break;
    case NORMALBUFFER:
        setColoursByNormalBuffer();
    default: // ISOPHOTES or REFLECTION
        break;
    }
}

int checkIfClosed(Vertex v) {
    HalfEdge *out = v.out;
    HalfEdge *currentEdge = out->twin->next;
    while(currentEdge != 0 && currentEdge != out) {
        currentEdge = currentEdge->twin->next;
    }
    return (currentEdge == out);
}

// returns vertex normals currently used in view
QVector<QVector3D> MainView::getNormals() {
    return vertexNormals;
}

void MainView::updateMeshBuffers() {
    setBuffers();
}

void MainView::updateMatrices() {

    modelViewMatrix.setToIdentity();
    modelViewMatrix.translate(QVector3D(0.0, 0.0, -3.0));
    modelViewMatrix.scale(scale);
    modelViewMatrix.rotate(rotAngleX, QVector3D(1.0, 0.0, 0.0));
    modelViewMatrix.rotate(rotAngleY, QVector3D(0.0, 1.0, 0.0));

    projectionMatrix.setToIdentity();
    projectionMatrix.ortho(-dispRatio, dispRatio, -1.0f, 1.0f, 0.1, 10.0);

//    float FoV;
//    float dispRatio;
//    float nearPlane;
//    float farPlane;

//    projectionMatrix.perspective(FoV, dispRatio, nearPlane, farPlane);

    normalMatrix = modelViewMatrix.normalMatrix();

    uniformUpdateRequired = true;
    update();

}

void MainView::updateUniforms() {

    // mainShaderProg should be bound at this point!

    glUniformMatrix4fv(uniModelViewMatrix, 1, false, modelViewMatrix.data());
    glUniformMatrix4fv(uniProjectionMatrix, 1, false, projectionMatrix.data());
    glUniformMatrix3fv(uniNormalMatrix, 1, false, normalMatrix.data());
    glUniform1f(uniIsoFreq, isophotesFrequency);

}

// Given a mesh, copy all properties
void MainView::setMesh(Mesh *mesh) {
    vertCount = mesh->Vertices.size();
    setCoords(mesh);
    setNormals(mesh);
    setPolyIndices(mesh);

    setLerpSlerpDiffs(mesh);
    setSurfaceSubdivDiffs(mesh);
    setLeftRightDiffs(mesh);
    setLimitDiffs(mesh);
    setVertsRegular(mesh);
    setBlendWeights(mesh);

    setColours();
}

void MainView::updateView() {
    addShaders();
    updateMeshBuffers();
    updateMatrices();
}

void MainView::setDrawMode(DrawMode dm) {
    drawMode = dm;
}

void MainView::setColourMode(ColourMode cm) {
    colourMode = cm;
    setColours();
    addShaders();
}


// ---

void MainView::initializeGL() {

    initializeOpenGLFunctions();
    qDebug() << ":: OpenGL initialized";

    debugLogger = new QOpenGLDebugLogger();
    connect( debugLogger, SIGNAL( messageLogged( QOpenGLDebugMessage ) ), this, SLOT( onMessageLogged( QOpenGLDebugMessage ) ), Qt::DirectConnection );

    if ( debugLogger->initialize() ) {
        qDebug() << ":: Logging initialized";
        debugLogger->startLogging( QOpenGLDebugLogger::SynchronousLogging );
        debugLogger->enableMessages();
    }

    QString glVersion;
    glVersion = reinterpret_cast<const char*>(glGetString(GL_VERSION));
    qDebug() << ":: Using OpenGL" << qPrintable(glVersion);

    // Enable depth buffer
    glEnable(GL_DEPTH_TEST);
    // Default is GL_LESS
    glDepthFunc(GL_LEQUAL);

    glEnable(GL_PRIMITIVE_RESTART);
    maxInt = ((unsigned int) -1);
    glPrimitiveRestartIndex(maxInt);

    // ---

    createShaderPrograms();

    addShaders();

    createBuffers();

    // ---

    updateMatrices();
}

void MainView::resizeGL(int newWidth, int newHeight) {

    qDebug() << ".. resizeGL";

    dispRatio = (float)newWidth/newHeight;
    updateMatrices();

}

void MainView::paintGL() {

    if (modelLoaded) {

        //glClearColor(0.0, 0.0, 0.0, 1.0);
        glClearColor(1.0, 1.0, 1.0, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        mainShaderProg->bind();

        if (uniformUpdateRequired) {
            updateUniforms();
            uniformUpdateRequired = false;
        }

        renderMesh();

        mainShaderProg->release();

    }
}

// ---

void MainView::renderMesh() {

    glBindVertexArray(meshVAO);

    enum DrawMode {WIREFRAME, FILLED, BOTH};

    if (drawMode != WIREFRAME) {
        glUniform1ui(uniShadingType, (drawMode == 2 ? 3 : 1) ); // 1 or 2, mapped to 1 or 3
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
        glDrawElements(GL_TRIANGLE_FAN, meshIBOSize2, GL_UNSIGNED_INT, 0);
    }

    if (drawMode != FILLED) {
        glUniform1ui(uniShadingType, drawMode); // 0 or 2
        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
        glDrawElements(GL_LINE_LOOP, meshIBOSize2, GL_UNSIGNED_INT, 0);
    }

    glBindVertexArray(0);
}

// Save current angle on mouse events
void MainView::setCurrentAngle(float angleX, float angleY) {
    currentAngleX = angleX;
    currentAngleY = angleY;
}

// Update rotate angle
void MainView::setRotateAngle(float angleX, float angleY) {
    rotAngleX = angleX;
    rotAngleY = angleY;
}

//// Update scale
//void MainView::setScale(double s) {
//    otherView->scale = s;
//    otherView->otherView->scale = s;

//    updateMatrices();
//    otherView->updateMatrices();
//    otherView->otherView->updateMatrices();
//}

// Mouse event catchers
void MainView::mousePressEvent(QMouseEvent* event) {
    setFocus();
    currentPos = event->pos();

    setCurrentAngle(rotAngleX, rotAngleY);
    otherView->setCurrentAngle(rotAngleX, rotAngleY);
    otherView->otherView->setCurrentAngle(rotAngleX, rotAngleY);
}

void MainView::mouseReleaseEvent(QMouseEvent* event) {
    currentPos = event->pos();

    setCurrentAngle(rotAngleX, rotAngleY);
    otherView->setCurrentAngle(rotAngleX, rotAngleY);
    otherView->otherView->setCurrentAngle(rotAngleX, rotAngleY);
}

void MainView::wheelEvent(QWheelEvent* event) {
    scale += event->delta() * 0.0001;

    otherView->scale = scale;
    otherView->otherView->scale = scale;

    updateMatrices();
    otherView->updateMatrices();
    otherView->otherView->updateMatrices();

}

void MainView::mouseMoveEvent(QMouseEvent* event) {
    QPointF diff = (event->pos() - currentPos) * 0.5;
    rotAngleX = fmod(currentAngleX + diff.y(), 360.0);
    rotAngleY = fmod(currentAngleY + diff.x(), 360.0);

    otherView->setRotateAngle(rotAngleX, rotAngleY);
    otherView->otherView->setRotateAngle(rotAngleX, rotAngleY);

    updateMatrices();
    otherView->updateMatrices();
    otherView->otherView->updateMatrices();
}

void MainView::onMessageLogged( QOpenGLDebugMessage Message ) {
    if (Message.id() == 131185 || Message.id() == 1285) {return;}
    qDebug() << " → Log:" << Message;
}
