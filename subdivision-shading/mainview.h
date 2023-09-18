#ifndef MAINVIEW_H
#define MAINVIEW_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions_4_1_Core>
#include <QOpenGLDebugLogger>

#include <QOpenGLShaderProgram>

#include <QMouseEvent>
#include "mesh.h"
#include "subdivider.h"

class MainView : public QOpenGLWidget, protected QOpenGLFunctions_4_1_Core {

    Q_OBJECT

public:
    MainView(QWidget *Parent = 0);
    ~MainView();

    // Several enums for draw-options
    enum DrawMode {WIREFRAME, FILLED, BOTH};
    enum NormalMode {FACEAVG, SUBDIV};
    enum ColourMode {REFLECTION, NORMALSDIFF, NORMALBUFFER, NORMALSDIFF_MESH, LIMIT_NORMALS_DIFF, EXTRORDVERTS, ISOPHOTES, LEFT_RIGHT_DIFF, BLENDWEIGHTS};

    // Subdivision scheme currently used
    Mesh::SubdivType subdivType;

    // Wireframe or Filled
    DrawMode drawMode;
    ColourMode colourMode;

    // Should surface and subdivision normals be blended based on blend weights
    bool blendNormals;

    bool modelLoaded;

    // The other view this view is compared to
    MainView *otherView;

    // Frequency of the isophotes
    int isophotesFrequency;

    // Several variables for scaling and rotation
    float scale;
    float dispRatio;
    float rotAngleY, rotAngleX, currentAngleX, currentAngleY;

    // Position of the mouse
    QPointF currentPos;

    // Should shader-uniforms be updated?
    bool uniformUpdateRequired;

    void updateMatrices();
    void updateUniforms();
    void updateMeshBuffers();
    void addShaders();

    // Set max difference value for scaling difference colors
    void setMaxDiff(float value);

    void updateView();

    // Given a mesh, copy all properties
    void setMesh(Mesh *mesh);

    void setDrawMode(DrawMode dm);
    void setColourMode(ColourMode cm);

    // Save current angle on mouse events
    void setCurrentAngle(float angleX, float angleY);

    // Update rotate angle
    void setRotateAngle(float angleX, float angleY);

    // Update scale
    void setScale(double s);

    // returns vertex normals currently used in view
    QVector<QVector3D> getNormals();
protected:
    void initializeGL();
    void resizeGL(int newWidth, int newHeight);
    void paintGL();

    unsigned int maxInt;

    void renderMesh();

    void mousePressEvent(QMouseEvent* event);
    void wheelEvent(QWheelEvent* event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
private:
    QOpenGLDebugLogger* debugLogger;

    QMatrix4x4 modelViewMatrix, projectionMatrix;
    QMatrix3x3 normalMatrix;

    // Uniforms
    GLint uniModelViewMatrix, uniProjectionMatrix, uniNormalMatrix, uniIsoFreq, uniShadingType;

    // ---

    QOpenGLShaderProgram* mainShaderProg;
    QOpenGLShader *vertReflectionShader, *reflectionFragShader;
    QOpenGLShader *vertshader, *fixedNormalFragShader, *isophotesFragShader;

    GLuint meshVAO, meshCoordsBO, meshNormalsBO, meshColoursBO, meshIndexBO;
    unsigned int meshIBOSize, meshIBOSize2;
    float maxDiff;
    int vertCount;

    // ---

    QVector<QVector3D> vertexCoords, vertexCoords2;
    QVector<QVector3D> vertexNormals, vertexNormals2;
    QVector<QVector3D> vertexColours, vertexColours2;
    QVector<float> lerpSlerpDiffs, faceSlerpDiffs, leftRightDiffs, limitDiffs, limitDiffs2, blendWeights, blendWeights2;
    QVector<bool> vertIsRegular;
    QVector<unsigned int> polyIndices, polyIndices2;

    void createShaderPrograms();
    void createBuffers();

    void setNormalsByFaces();
    void setBuffers();
    void setNormals(Mesh *mesh);
    void calcTestNormals(Mesh* mesh);

    void setDefaultColors();

    // Different colour modes. Set colour of vertex based on property
    void setColoursByFaceNormals();
    void setColoursBySubdividedNormals();
    void setColoursByNormalBuffer();
    void setColoursByDiffLerpSlerp();
    void setColoursByDiffSurfaceSubdiv();
    void setColoursByDiffLimit();
    void setColoursByExtrordVerts();
    void setColoursByDiffLeftRight();
    void setColoursByBlendWeights();
    void setColours2(Mesh *mesh);

    // Copy coords and other properties from given mesh
    void setCoords(Mesh *mesh);
    void setPolyIndices(Mesh *mesh);
    void setColours();
    void setLimitDiffs(Mesh *mesh);
    void setLeftRightDiffs(Mesh *mesh);
    void setSurfaceSubdivDiffs(Mesh *mesh);
    void setProperties(Mesh *mesh);

    // Set differences between LinEar inteRPolation (LERP) and Spherical Linear inteRPolation (SLERP) of normals
    void setLerpSlerpDiffs(Mesh *mesh);

    // Set property of vertices: regular or not
    void setVertsRegular(Mesh *mesh);

    // Set the blend weights of the vertices
    void setBlendWeights(Mesh *mesh);
private slots:
    void onMessageLogged( QOpenGLDebugMessage Message );
};

#endif // MAINVIEW_H
