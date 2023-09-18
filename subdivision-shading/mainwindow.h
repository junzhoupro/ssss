#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "objfile.h"
#include <QFileDialog>
#include "mesh.h"
#include "loopsubdivider.h"
#include "catclarksubdivider.h"
#include "mainview.h"
#include <QSplitter>
#include <QtWidgets/QLabel>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {

    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    // Subdivision scheme currently used
    Mesh::SubdivType subdivType;

    // Option for normal blending: no blending, linear blending, Subdivision Blending
    enum BlendOption {NOBLEND, LINBLEND, SUBDIVBLEND, REALNOBLEND};
    enum SharpOption {REGULAR, BOUNDARY, ZERO};

    BlendOption blendOption;

    // Meshes for Loop subdivision
    QVector<Mesh> loopMeshes;

    // Meshes for Catmull-Clark subdivision
    QVector<Mesh> catclarkMeshes;

    // Meshes currently used
    QVector<Mesh> *Meshes;

    // Function for loading a new object file
    void loadOBJ();
    void loadSharpness();
    void updateModel(OBJFile newModel, bool reset);


private slots:
    // Event catchers
    void on_LoadOBJ_clicked();

    void on_LoadSharpness_clicked();

    void on_wireframeOption_clicked();

    void on_filledOption_clicked();

    void on_fixedColour_clicked();

    void on_normalBuffer_clicked();

    void on_geomSubdiv_valueChanged(int arg1);

    void on_slerpIterations_valueChanged(int arg1);

    void on_loop_clicked();

    void on_cc_clicked();

    void on_modButterOption_clicked();

    void on_normalSubdivStart_left_valueChanged(int value);

    void on_normalSubdivStop_left_valueChanged(int value);

    void on_normalSubdivision_left_toggled(bool checked);

    void on_coordsToLimit_toggled(bool checked);

    void on_normalsToLimit_toggled(bool checked);

    void on_diffLimitNormals_clicked();

    void on_normalSubdivision_right_toggled(bool checked);

    void on_normalSubdivStart_right_valueChanged(int arg1);

    void on_normalSubdivStop_right_valueChanged(int arg1);

    void on_slerpIterations_right_valueChanged(int arg1);

    void on_extraordinaryVerts_pressed();

    void on_isophotes_clicked();

    void on_isoFreq_valueChanged(int arg1);

    void on_blendWeights_clicked();

    void on_blendNormalsFrom_valueChanged(int arg1);

    void on_subdivBlending_clicked();

    void on_linBlending_clicked();

    void on_noBlending_clicked();

    void on_p_valueChanged(int value);

    void on_radioButton_toggled(bool checked);

    void on_radioButton_2_toggled(bool checked);

    void on_radioButton_3_toggled(bool checked);

    void on_radioButton_4_toggled(bool checked);

    void on_bothOption_clicked();

    void on_sharpDarts_toggled(bool checked);

    void on_sharpRegular_clicked();

    void on_sharpBoundary_clicked();

    void on_sharpZero_clicked();

    void on_screenShot_clicked();

//    void on_doubleSpinBox_valueChanged(int arg1);


private:
    Ui::MainWindow *ui;

    // LoopSubdivider or CatClarkSubdivider
    LoopSubdivider loopDivider;
    CatClarkSubdivider catclarkDivider;
    Subdivider *divider;

    // Normal subdivision enabled for left/right view?
    // Should coordinates be projected to limit after subdivision?
    // Should normals be projected to limit after subdivision?
    bool normalSubdiv_left, normalSubdiv_right, coordsToLimit, normalsToLimit;

    // Number of geometry subdivision steps
    // At what Mesh level should blend values be initialized?
    int geomSubdivSteps, blendLevel;
    float p;

    // There are three different subdivision blending approaches (orig = setting EVs to 1, ours A = setting EVs to 1/lim, ours B = init one-ring NBH to 1)
    short subdivApproach = 2;

    // At what mesh index should normals be subdivided? 1 for immediate subdivision (based on index 0)
    // At what mesh index should normal subdivision stop (should be >= geomSubdivSteps)
    // Number of SLERP iterations
    int normalSubdivStart_left, normalSubdivStop_left, slerpIterations;
    int normalSubdivStart_right, normalSubdivStop_right, slerpIterations_2;

    // maximum difference in normals between left and right view (for color scaling)
    float maxDiff;

    // reset ui for loading new model
    void resetControls();

    void switchSubdivisionType(Mesh::SubdivType type);

    // New data for view: update buffers
    void updateLeftView();
    void updateRightView();
    void updateDiffView();

    Mesh *getMesh(short level);

    // Subdivide at least 'level' times
    void ensureSubdivision(int level);

    float getSlerpAngle(int normalStop, int geomSteps, int slerpIt);

    // NORMALBUFFER, NORMALSDIFF_MESH, LIMIT_NORMALS_DIFF, EXTRORDVERTS, ISOPHOTES, LEFT_RIGHT_DIFF, BLENDWEIGHTS
    void setColourMode(MainView::ColourMode colourMode);
    void setIsoFrequency(int arg1);
    void updateViews(bool left, bool right);
    void subdivBlend();
    void showSlerpAngle(short normalSDstop, short slerpIt, QLabel *label);
    void updateView(bool normalSD, short normalSDstart, short normalSDstop, short slerpIt, MainView *view, QLabel *slerpAngleLabel);
    void setMaxDiff(float value);

    int sharpEV = SharpOption::REGULAR;
    bool sharpEnd = false;
    bool noSharpness = false;
    OBJFile currentModel;
};

#endif // MAINWINDOW_H
