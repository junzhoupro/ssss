#include "objfile.h"

#include <QDebug>
#include <QFile>


OBJFile::OBJFile() {

}

OBJFile::OBJFile(QString fileName) {
    qDebug() << "✓✓ OBJFile constructor";

    qDebug() << ":: Loading" << fileName;
    QFile newModel(fileName);

    fileLoaded = newModel.open(QIODevice::ReadOnly);

    if(fileLoaded) {
        QTextStream fileContents(&newModel);

        QString currentLine;
        QStringList values;
        QStringList indices;

        unsigned short k;

        vertexCoords.clear();
        textureCoords.clear();
        vertexNormals.clear();
        faceCoordInd.clear();
        faceTexInd.clear();
        faceNormalInd.clear();

        int offset = 1;

        while(!fileContents.atEnd()) {

            currentLine = fileContents.readLine();
            values = currentLine.split(" ");


            if(values[0] == "OFF" && values[1] == "0") {
                offset = 0;
            }
            else if (values[0] == "v") {
                // qDebug() << "Vertex coords";
                // Only x, y and z. If there's a w value (homogenous coordinates), ignore it.
                vertexCoords.append(QVector3D(values[1].toFloat(), values[2].toFloat(), values[3].toFloat() ));
            }
            else if (values[0] == "vt") {
                // qDebug() << "Texture coords";
                // Only u and v. If there's a w value (barycentric coordinates), ignore it, it can be retrieved from 1-u-v.
                textureCoords.append(QVector2D(values[1].toFloat(), values[2].toFloat() ));
            }
            else if (values[0] == "vn") {
                // qDebug() << "Vertex normal";
                vertexNormals.append(QVector3D(values[1].toFloat(), values[2].toFloat(), values[3].toFloat() ));
            }
            else if (values[0] == "fs") {
                QVector<unsigned int> fs;
                for (k=1; k<values.size(); k++) {
                    fs.append(values[k].toInt());
                }
                faceSharpness.append(fs);
            }
            else if (values[0] == "as") {
                for(int i = 0; i < values[1].toInt(); i++) {
                    QVector<unsigned int> fs;
                    for (k=2; k<values.size(); k++) {
                        fs.append(values[k].toInt());
                    }
                    faceSharpness.append(fs);
                }
            }
            else if (values[0] == "f") {
                // qDebug() << "Face";

                for (k=1; k<values.size(); k++) {
                    indices = values[k].split("/");

                    // Note -1, OBJ starts indexing from 1.

                    faceCoordInd.append(indices[0].toInt() - offset );

                    if (indices.size() > 1) {
                        if (!indices[1].isEmpty()) {
                            faceTexInd.append(indices[1].toInt() - offset );
                        }

                        if (indices.size() > 2) {
                            if (!indices[2].isEmpty()) {
                                faceNormalInd.append(indices[2].toInt() - offset );
                            }
                        }
                    }

                }

                faceValences.append(k-1);

            }
            else {
                qDebug() << " * Line contents ignored," << currentLine;
            }

        }

        newModel.close();

    }

}

bool OBJFile::loadedFile() {
    return fileLoaded;
}

OBJFile::~OBJFile() {
    qDebug() << "✗✗ OBJFile destructor";
}

void OBJFile::updateSharpness(QVector<QVector<unsigned int>> newSharpness) {
    qDebug() << faceNormalInd;
    qDebug() << faceSharpness;
    for(int i = 0; i < faceSharpness.size(); i++) {
        qDebug() << i;
        faceSharpness[i].clear();
        faceSharpness[i].squeeze();
    }
    faceSharpness.clear();
    faceSharpness.squeeze();
    faceSharpness = newSharpness;
}
