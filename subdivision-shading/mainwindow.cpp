#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "assert.h"
#include <math.h>
#include <QDateTime>


int max(int a, int b) {
    return (b > a ? b : a);
}

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow) {
    qDebug() << "✓✓ MainWindow constructor";
    ui->setupUi(this);

     //Default subdivision scheme is Loop
        divider = &(loopDivider);
        Meshes = &loopMeshes;
        subdivType = Mesh::SubdivType::LOOP;
//    divider = &(catclarkDivider);
//    Meshes = &catclarkMeshes;
//    subdivType = Mesh::SubdivType::CATCLARK;

    // init values to defaults
    geomSubdivSteps = 0;
    normalSubdivStart_left = 1;
    normalSubdivStop_left = 1;
    slerpIterations = 0;
    normalSubdiv_left = false;
    coordsToLimit = false;
    normalsToLimit = false;
    blendLevel = 2;

    normalSubdivStart_right = 1;
    normalSubdivStop_right = 1;
    slerpIterations_2 = 0;
    normalSubdiv_right = false;
    maxDiff = 1.0;

    // Couple views
    ui->LeftDisplay->otherView = ui->RightDisplay;
    ui->RightDisplay->otherView = ui->DiffDisplay;
    ui->DiffDisplay->otherView = ui->LeftDisplay;
}

MainWindow::~MainWindow() {
    qDebug() << "✗✗ MainWindow destructor";
    delete ui;

    loopMeshes.clear();
    loopMeshes.squeeze();

    catclarkMeshes.clear();
    catclarkMeshes.squeeze();
}

void setGroupEnabled(QGroupBox *box) {
    QObjectList children = box->children();
    int cnt = children.length();
    box->setEnabled(true);
    for(int i = 0; i < cnt; i++) {
        QRadioButton *button = (QRadioButton *) children[i];
        button->setEnabled(true);
    }
}

void MainWindow::resetControls() {
    // reset variables
    geomSubdivSteps = 0;
    normalSubdivStart_left = 1;
    normalSubdivStart_right = 1;
    normalSubdivStop_left = 1;
    normalSubdivStop_right = 1;
    blendLevel = 0;
    p = 1.0;
    blendOption = BlendOption::REALNOBLEND;

    // reset UI
    ui->geomSubdiv->blockSignals(true);
    ui->geomSubdiv->blockSignals(true);
    ui->geomSubdiv->blockSignals(true);
    ui->geomSubdiv->blockSignals(true);
    ui->geomSubdiv->blockSignals(true);
    ui->linBlending->blockSignals(true);
    ui->subdivBlending->blockSignals(true);
    ui->noBlending->blockSignals(true);

    ui->geomSubdiv->setValue(0);
    ui->normalSubdivStart_left->setValue(0);
    ui->normalSubdivStop_left->setValue(0);
    ui->normalSubdivStart_right->setValue(0);
    ui->normalSubdivStop_right->setValue(0);
    ui->linBlending->setEnabled(false);
    ui->subdivBlending->setEnabled(false);
    ui->noBlending->setEnabled(true);

    ui->geomSubdiv->blockSignals(false);
    ui->geomSubdiv->blockSignals(false);
    ui->geomSubdiv->blockSignals(false);
    ui->geomSubdiv->blockSignals(false);
    ui->geomSubdiv->blockSignals(false);
    ui->linBlending->blockSignals(false);
    ui->subdivBlending->blockSignals(false);
    ui->noBlending->blockSignals(false);

    // update both views
    updateViews(true, true);
}

// Updates the difference-view (difference between left and right)
void MainWindow::updateDiffView() {
    Mesh *m = getMesh(geomSubdivSteps);

    QVector<QVector3D> normalsLeft = ui->LeftDisplay->getNormals();
    QVector<QVector3D> normalsRight = ui->RightDisplay->getNormals();
    assert(normalsLeft.size() == normalsRight.size());

    // calculate maximum difference in normals between two views
    maxDiff = 0.0;
    for(int i = 0; i < normalsLeft.size(); i++) {
        m->Vertices[i].leftNormal = normalsLeft[i];
        m->Vertices[i].rightNormal = normalsRight[i];
        maxDiff = fmax(maxDiff, 1.0 - QVector3D::dotProduct(normalsLeft[i], normalsRight[i]));
    }

    // use same mesh (geometry)
    ui->DiffDisplay->setMesh(m);
    ui->DiffDisplay->setColourMode(MainView::ColourMode::LEFT_RIGHT_DIFF);
    ui->DiffDisplay->updateView();
}

// Sets colour modes for left and right view
void MainWindow::setColourMode(MainView::ColourMode colourMode) {
    ui->LeftDisplay->setColourMode( colourMode);
    ui->RightDisplay->setColourMode( colourMode);
    ui->LeftDisplay->updateView();
    ui->RightDisplay->updateView();
}

// Sets frequency of isophotes
void MainWindow::setIsoFrequency(int arg1) {
    ui->LeftDisplay->isophotesFrequency = arg1;
    ui->RightDisplay->isophotesFrequency = arg1;
    ui->LeftDisplay->updateView();
    ui->RightDisplay->updateView();
}

// Switches subdivision type if not the same as current
void MainWindow::switchSubdivisionType(Mesh::SubdivType type) {
    if(type != subdivType) {
        subdivType = type;
        if(type == Mesh::SubdivType::LOOP) {
            divider = &(loopDivider);
            Meshes = &loopMeshes;
        } else if(type == Mesh::SubdivType::CATCLARK) {
            divider = &(catclarkDivider);
            Meshes = &catclarkMeshes;
        }

        ui->LeftDisplay->subdivType = type;
        ui->RightDisplay->subdivType = type;

        // If Model is loaded
        if(Meshes->size()) {
            updateViews(true, true);
        }
    }
}

// Gets maximum angle between SLERP and LERP normals in radians and converts it to degrees
float MainWindow::getSlerpAngle(int normalStop, int geomSteps, int slerpIt) {
    float angle, deg;

    if(normalStop < geomSteps || slerpIt == 0) {
        deg = 0;
    } else {
        angle = getMesh(normalStop)->calcMaxSlerpAngle();
        qDebug() << "slerp angle: " << angle;
        deg = angle * 360.0 / (2.0 * M_PI);
//        deg = roundf(deg * 100.0) / 100.0;
    }

    qDebug() << "slerp angle in deg: " << deg;

    return deg;
}

// Update maximum angle between SLERP and LERP normal
void MainWindow::showSlerpAngle(short normalSDstop, short slerpIt, QLabel *label) {
//    char str[6];
//    sprintf(str, "%02.6f", getSlerpAngle(normalSDstop, geomSubdivSteps, slerpIt));
//    label->setText(str);
}

// Loads a new object picked by user via file dialog
void MainWindow::loadOBJ() {
    OBJFile newModel = OBJFile(QFileDialog::getOpenFileName(this, "Import OBJ File", "./../subdivision-shading_loop2/models", tr("Obj Files (*.obj)"))); // ../comparison/models/

    if(newModel.loadedFile()) {
        updateModel(newModel, true);
    }
}

void MainWindow::updateModel(OBJFile newModel, bool reset) {
    currentModel = newModel;
    // Clear meshes
    loopMeshes.clear();
    loopMeshes.squeeze();
    catclarkMeshes.clear();
    catclarkMeshes.squeeze();

    // Load new model
    loopMeshes.append( Mesh(&newModel) );
    catclarkMeshes.append( Mesh(&newModel) );

    // Sets surface normals
    loopMeshes[0].setFaceNormals();
    catclarkMeshes[0].setFaceNormals();

    Mesh *origMesh = getMesh(0);
    ui->LeftDisplay->setMesh(origMesh);
    ui->RightDisplay->setMesh(origMesh);
    ui->DiffDisplay->setMesh(origMesh);

    ui->LeftDisplay->modelLoaded = true;
    ui->RightDisplay->modelLoaded = true;
    ui->DiffDisplay->modelLoaded = true;

    if (reset) {
        resetControls();
    }

    setGroupEnabled(ui->subdivBox);
    setGroupEnabled(ui->subdivBox_left);
    setGroupEnabled(ui->subdivBox_right);
    setGroupEnabled(ui->drawModeBox);
    ui->normalSubdivision_left->setChecked(false);
    ui->normalSubdivision_right->setChecked(false);
    ui->normalSubdivStart_left->setEnabled(false);
    ui->normalSubdivStart_right->setEnabled(false);
    ui->normalSubdivStop_left->setEnabled(false);
    ui->normalSubdivStop_right->setEnabled(false);
    setGroupEnabled(ui->colourModeBox);
    setGroupEnabled(ui->blendNormalBox);

    updateViews(true, true);
}

// Loads a new sharpness preset picked by user via file dialog
void MainWindow::loadSharpness() {
    qDebug() << "Loading new sharpness preset";
    OBJFile newSharpness = OBJFile(QFileDialog::getOpenFileName(this, "Import Sharpness File", "/home/pieter/Dropbox/OpenGL_Code/Models", tr("Sharpness Files (*.sharp)"))); // ../comparison/models/

    if(newSharpness.loadedFile()) {
        currentModel.updateSharpness(newSharpness.faceSharpness);
        updateModel(currentModel, false);
    }
}

Mesh *MainWindow::getMesh(short level) {
    return &((*Meshes)[level]);
}

void MainWindow::on_LoadOBJ_clicked() {
    loadOBJ();
}

void MainWindow::on_LoadSharpness_clicked() {
    loadSharpness();
}

// based on new settings (for left or right view), subdivide and compute normals
void MainWindow::updateViews(bool left, bool right) {
    // make sure the mesh is at least subdivided 'geomSubdivSteps' times
    ensureSubdivision(geomSubdivSteps);

    if(coordsToLimit) {
        divider->calculateLimitCoords(getMesh(geomSubdivSteps));
        getMesh(geomSubdivSteps)->setLimitCoords();
    } else {
        getMesh(geomSubdivSteps)->setOrigCoords();
    }

    // Blending of normals (linear or Subdivision Blending)
    if(blendOption == BlendOption::REALNOBLEND) {qDebug() << "blendOption == BlendOption::REALNOBLEND";
        for(short i = 0; i <= geomSubdivSteps; i++) {
            divider->calculateNoBlendWeights(getMesh(i));
        }
    } else if(blendOption == BlendOption::LINBLEND) {qDebug() << "blendOption == BlendOption::LINBLEND";
        divider->setLinearBlendWeights(getMesh(blendLevel), sharpEV, sharpEnd);
        for(short i = blendLevel+1; i <= geomSubdivSteps; i++) {
            divider->calculateLinearBlendWeights(getMesh(i));
        }
    } else if(blendOption == BlendOption::SUBDIVBLEND) {qDebug() << "blendOption == BlendOption::SUBDIVBLEND" << "sharpEV" << sharpEV;
        subdivBlend();
    }

    if(left) { // If settings for left view changed
        updateView(normalSubdiv_left, normalSubdivStart_left, normalSubdivStop_left, slerpIterations, ui->LeftDisplay, ui->slerpAngleLabel);
    }
    if(right) { // If settings for left view changed
        updateView(normalSubdiv_right, normalSubdivStart_right, normalSubdivStop_right, slerpIterations_2, ui->RightDisplay, ui->slerpAngleLabel_2);
    }
    if(left || right) { // if one of the views changed -> update difference view
        //updateDiffView();
    }
}

void MainWindow::subdivBlend() {

  divider->setInitialBlendWeights(getMesh(blendLevel), true, subdivApproach, sharpEV, sharpEnd);
  for(short i = blendLevel+1; i <= geomSubdivSteps; i++) {
      divider->calculateBlendWeights(getMesh(i), subdivApproach, sharpEV, sharpEnd);
  }
  qDebug() << "p=" << p;
  divider->calculateLimitBlendWeights(getMesh(geomSubdivSteps), p);
}

// Update one of the views (left or right)
void MainWindow::updateView(bool normalSD, short normalSDstart, short normalSDstop, short slerpIt, MainView *view, QLabel *slerpAngleLabel) {
    Mesh *mesh, *geomMesh;
    int faceLevel;

    // get right mesh for geometry
    geomMesh = getMesh(geomSubdivSteps);

    // determine level where surface normals should be computed
    if(normalSDstop < geomSubdivSteps || !normalSD) {
        // if not normal subdivision -> surface normals at last level
        faceLevel = geomSubdivSteps;
    } else {
        faceLevel = normalSDstart-1;
    }

    // set surface normals at right level
    // v->normal = v->surfaceNormal;
    // v->sharpNormals = v->sharpSurfaceNormals;
    getMesh(faceLevel)->setFaceNormals();

    // if normal subdivision
    if(normalSD && normalSDstop >= geomSubdivSteps) {
        // make sure that the mesh is subdivided enough (for topology of normal subdivision)
        ensureSubdivision(normalSDstop); qDebug() << "######";
        for(short k = normalSDstart; k <= normalSDstop; k++) {qDebug() << "k" << k << subdivApproach << sharpEV;
            mesh = getMesh(k);
//            mesh->substractEV(0);
//            qDebug() << "subdivType" << subdivType;
            divider->calculateNormals(mesh, slerpIt, blendOption, subdivApproach, sharpEV, subdivType);
            // v->normal = v->sphereNormal;
            // v->sharpNormals = v->sharpSphereNormals;
            mesh->setSubdivideNormals();

        }

        if(normalSDstop > geomSubdivSteps) {
            geomMesh->setNormalsAtDepth(normalSDstop-geomSubdivSteps);
        }

        if(blendOption != BlendOption::REALNOBLEND) {qDebug() << "BLEND!";
//        if(blendOption != BlendOption::NOBLEND) {
            // v->normal = ;
            // v->sharpNormals = ;
            geomMesh->blendNormals();
        }
    }

    if(normalsToLimit) { // project normals to limit
        divider->calculateLimitNormals(geomMesh);
        geomMesh->setLimitNormals();
    }

    // update buffers
    view->setMesh(geomMesh);
    view->updateView();

    // show maximum angle between LERP and SLERP normals
    showSlerpAngle(normalSDstop, slerpIt, slerpAngleLabel);

}

// Make sure the mesh is subdivided at least 'level' times
void MainWindow::ensureSubdivision(int level) {
    for(short k = Meshes->size(); k <= level; k++) {
        Meshes->append(Mesh());
        divider->subdivide(getMesh(k-1), getMesh(k));
        divider->calculateCoords(getMesh(k));
        qDebug() << "vertice size" << getMesh(k)->Vertices.size();
    }
}

// Event catchers
void MainWindow::on_loop_clicked() {
    switchSubdivisionType(Mesh::SubdivType::LOOP);
}

void MainWindow::on_cc_clicked() {
    switchSubdivisionType(Mesh::SubdivType::CATCLARK);
}

void MainWindow::on_modButterOption_clicked() {
    switchSubdivisionType(Mesh::SubdivType::MODBUTTERFLY);
}

void MainWindow::setMaxDiff(float value) {
    ui->LeftDisplay->setMaxDiff(value*maxDiff);
    ui->RightDisplay->setMaxDiff(value*maxDiff);
    ui->DiffDisplay->setMaxDiff(value*maxDiff);
    updateViews(true, true);
}

// DRAW MODES
void MainWindow::on_wireframeOption_clicked() {
    ui->LeftDisplay->setDrawMode(MainView::DrawMode::WIREFRAME);
    ui->RightDisplay->setDrawMode(MainView::DrawMode::WIREFRAME);
    ui->DiffDisplay->setDrawMode(MainView::DrawMode::WIREFRAME);

    updateViews(true, true);
}

void MainWindow::on_filledOption_clicked() {
    ui->LeftDisplay->setDrawMode(MainView::DrawMode::FILLED);
    ui->RightDisplay->setDrawMode(MainView::DrawMode::FILLED);
    ui->DiffDisplay->setDrawMode(MainView::DrawMode::FILLED);

    updateViews(true, true);
}

// COLOUR MODES
void MainWindow::on_fixedColour_clicked() {
    setColourMode( MainView::ColourMode::REFLECTION);
}

void MainWindow::on_diffLimitNormals_clicked() {
    setColourMode( MainView::ColourMode::LIMIT_NORMALS_DIFF);
}

void MainWindow::on_normalBuffer_clicked() {
    setColourMode( MainView::ColourMode::NORMALBUFFER);;
}

void MainWindow::on_extraordinaryVerts_pressed() {
    setColourMode(MainView::ColourMode::EXTRORDVERTS);
}

void MainWindow::on_isophotes_clicked() {
    setColourMode( MainView::ColourMode::ISOPHOTES);
}

void MainWindow::on_isoFreq_valueChanged(int arg1) {
    setIsoFrequency(arg1);
}

void MainWindow::on_blendWeights_clicked() {
    setColourMode( MainView::ColourMode::BLENDWEIGHTS);
}


// SUBDIVISION
void MainWindow::on_geomSubdiv_valueChanged(int value) {
    geomSubdivSteps = value;
    ui->DiffDisplay->setDrawMode(MainView::DrawMode::FILLED);
    updateViews(true, true);
}

void MainWindow::on_slerpIterations_valueChanged(int value) {
    slerpIterations = value;

    if(normalSubdiv_left) {
        updateViews(true, false);
    }
}

void MainWindow::on_slerpIterations_right_valueChanged(int arg1) {
    slerpIterations_2 = arg1;

    if(normalSubdiv_right) {
        updateViews(false, true);
    }
}

void MainWindow::on_normalSubdivStart_left_valueChanged(int value) {
    normalSubdivStart_left = value;
    if(normalSubdivStart_left > normalSubdivStop_left) {
        ui->normalSubdivStop_left->setValue(value);
    }
    updateViews(true, false);
}

void MainWindow::on_normalSubdivStart_right_valueChanged(int arg1) {
    normalSubdivStart_right = arg1;
    if(normalSubdivStart_right > normalSubdivStop_right) {
        ui->normalSubdivStop_right->setValue(arg1);
    }
    updateViews(false, true);
}

void MainWindow::on_normalSubdivStop_left_valueChanged(int value) {
    normalSubdivStop_left = value;
    if(normalSubdivStart_left > normalSubdivStop_left) {
        ui->normalSubdivStart_left->setValue(value);
    }
    updateViews(true, false);
}

void MainWindow::on_normalSubdivStop_right_valueChanged(int arg1) {
    normalSubdivStop_right = arg1;
    if(normalSubdivStart_right > normalSubdivStop_right) {
        ui->normalSubdivStart_right->setValue(arg1);
    }
    updateViews(false, true);
}

void MainWindow::on_normalSubdivision_left_toggled(bool checked) {
    normalSubdiv_left = checked;
    ui->normalSubdivStart_left->setEnabled(checked);
    ui->normalSubdivStop_left->setEnabled(checked);
    updateViews(true, false);
}


void MainWindow::on_normalSubdivision_right_toggled(bool checked) {
    normalSubdiv_right = checked;subdivBlend();
    ui->normalSubdivStart_right->setEnabled(checked);
    ui->normalSubdivStop_right->setEnabled(checked);
    updateViews(false, true);
}

void MainWindow::on_coordsToLimit_toggled(bool checked) {
    coordsToLimit = checked;
    updateViews(true, true);
}

void MainWindow::on_blendNormalsFrom_valueChanged(int arg1) {
    blendLevel = arg1;
    updateViews(true, true);
}

void MainWindow::on_normalsToLimit_toggled(bool checked) {
    normalsToLimit = checked;
    updateViews(true, true);
}

void MainWindow::on_subdivBlending_clicked() {
    blendOption = SUBDIVBLEND;
    blendLevel = ui->blendNormalsFrom->value();
    ui->LeftDisplay->blendNormals = true;
    ui->RightDisplay->blendNormals = true;
    updateViews(true, true);
}

void MainWindow::on_linBlending_clicked() {
    blendOption = LINBLEND;
    blendLevel = ui->blendNormalsFrom->value();
    ui->LeftDisplay->blendNormals = true;
    ui->RightDisplay->blendNormals = true;
    updateViews(true, true);
}

void MainWindow::on_noBlending_clicked() {
    blendOption = BlendOption::REALNOBLEND;
    blendLevel = ui->blendNormalsFrom->value();
    ui->LeftDisplay->blendNormals = false;
    ui->RightDisplay->blendNormals = false;
    updateViews(true, true);
}

void MainWindow::on_sharpRegular_clicked() {
    sharpEV = SharpOption::REGULAR;
    updateViews(true, true);
}

void MainWindow::on_sharpBoundary_clicked() {
    sharpEV = SharpOption::BOUNDARY;
    updateViews(true, true);
}

void MainWindow::on_sharpZero_clicked() {
    sharpEV = SharpOption::ZERO;
    updateViews(true, true);
}

void MainWindow::on_p_valueChanged(int tentimesP) {
    p = ((double) tentimesP) / 10.0;
    char str[10];
    sprintf(str, "p : %01.1f", p);
    ui->pLabel->setText(QString::fromStdString(str));
    updateViews(true, true);
    //updateViews(true, true); // Why twice?
}

void MainWindow::on_radioButton_toggled(bool checked) {
  subdivApproach = 0;
  updateViews(true, true);
  //subdivBlend();
}

void MainWindow::on_radioButton_2_toggled(bool checked) {
  subdivApproach = 1;
  updateViews(true, true);
}

void MainWindow::on_radioButton_3_toggled(bool checked) {
  subdivApproach = 2;
  updateViews(true, true);
}

void MainWindow::on_radioButton_4_toggled(bool checked) {
  subdivApproach = 3;
  updateViews(true, true);
}

void MainWindow::on_bothOption_clicked() {
  ui->LeftDisplay->setDrawMode(MainView::DrawMode::BOTH);
  ui->RightDisplay->setDrawMode(MainView::DrawMode::BOTH);
  ui->DiffDisplay->setDrawMode(MainView::DrawMode::BOTH);

  updateViews(true, true);
}

void MainWindow::on_sharpDarts_toggled(bool checked) {
    sharpEnd = checked;
    updateViews(true, true);
}

void MainWindow::on_screenShot_clicked() {
    // a function that returns the current widget.
    QWidget* widget = ui->LeftDisplay;
    if (not widget)
        return;

    QPixmap screenshot = widget->grab(widget->rect());
    if (screenshot.isNull())
        printf("Error printing widget");

    QDateTime d = QDateTime::currentDateTime();
    QString date = d.toString("dd.MM.yyyy.hh:mm:ss.zzz") + ".png";
    QString fileName =  "./../subdivision-shading_loop2/screenshots/" + date;

    bool s = screenshot.save(fileName, "PNG", 100);
    if (not s)
        qDebug() << "Error printing widget " << date;
    else
        qDebug() << "Sreenshot saved as " << date;
}


