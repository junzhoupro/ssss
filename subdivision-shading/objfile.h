#ifndef OBJFILE_H
#define OBJFILE_H

#include <QString>
#include <QVector>
#include <QVector2D>
#include <QVector3D>

class OBJFile {

public:
    OBJFile(QString fileName);
    OBJFile();
    ~OBJFile();

    QVector<QVector3D> vertexCoords;
    QVector<QVector2D> textureCoords;
    QVector<QVector3D> vertexNormals;
    QVector<unsigned short> faceValences;
    QVector<unsigned int> faceCoordInd;
    QVector<unsigned int> faceTexInd;
    QVector<unsigned int> faceNormalInd;
    QVector<QVector<unsigned int>> faceSharpness;

    void updateSharpness(QVector<QVector<unsigned int>> newSharpness);

    bool loadedFile();
private:
    bool fileLoaded;
};

#endif // OBJFILE_H
