#-------------------------------------------------
#
# Project created by QtCreator 2017-10-18T14:23:13
#
#-------------------------------------------------
LIBS    += -lSDL2 -lGLEW -lGL -lGLU
#-lglu32 -lopengl32 -lfreeglut
#-lGL -lglut -lGLEW -lglew32 -lopengl32 -lWs2_32 -lole32 -lcomctl32 -lgdi32 -lcomdlg32 -luuid

QT       = core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QMAKE_CXXFLAGS += -std=c++11

TARGET = 3DScan_Qt
TEMPLATE = app


#SOURCES += $$files(source/*.cpp)

SOURCES += $$files(source/Point3d.cpp)
SOURCES += $$files(source/KDTree.cpp)
SOURCES += $$files(source/GLcamera.cpp)
SOURCES += $$files(source/GLwidget.cpp)
SOURCES += $$files(source/MainWindow.cpp)
SOURCES += $$files(source/OpenGL_App_QT.cpp)


#HEADERS  += $$files(include/*.h)
HEADERS  += $$files(include/GLcamera.h)
HEADERS  += $$files(include/GLwidget.h)
HEADERS  += $$files(include/KDTree.h)
HEADERS  += $$files(include/MainWindow.h)
HEADERS  += $$files(include/Point3d.h)
