#-------------------------------------------------
#
# Project created by QtCreator 2019-10-09T08:47:57
#
#-------------------------------------------------

QT       += core gui printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = plasmaSolver
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

SOURCES += \
        crosssection.cpp \
        main.cpp \
        mainwindow.cpp \
        qcustomplot.cpp \
        simulationdata.cpp \
        simulationsolver.cpp \
        simulationtools.cpp \
        reaction.cpp

HEADERS += \
        crosssection.h \
        mainwindow.h \
        qcustomplot.h \
        simulationdata.h \
        simulationsolver.h \
        simulationtools.h \
        reaction.h

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

DISTFILES += \
    reactions/e+Ar_2e+Ar+.txt \
    reactions/e+Ar_e+Ar.txt \
    reactions/e+Ar_e+Ars.txt \
    reactions/e+Ars_2e+Ar+.txt \
    comsol/e+Ar_2e+Ar+.txt

OTHER_FILES += \
    comsol/e+Ar_e+Ar.txt \
    comsol/e+Ar_e+Ars.txt \
    comsol/e+Ars_2e+Ar+.txt
