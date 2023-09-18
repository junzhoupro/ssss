#-------------------------------------------------
#
# Project created by QtCreator 2016-12-01T23:50:32
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4.8): QT += widgets

TARGET = LoopSubdiv
TEMPLATE = app

SOURCES += main.cpp\
        mainwindow.cpp \
    mainview.cpp \
    objfile.cpp \
    mesh.cpp \
    subdivider.cpp \
    loopsubdivider.cpp \
    catclarksubdivider.cpp

HEADERS  += mainwindow.h \
    mainview.h \
    face.h \
    mesh.h \
    vertex.h \
    halfedge.h \
    objfile.h \
    subdivider.h \
    loopsubdivider.h \
    catclarksubdivider.h

FORMS    += mainwindow.ui

RESOURCES += \
    resources.qrc

#todo: limit normals
