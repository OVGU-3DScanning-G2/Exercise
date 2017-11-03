// OpenGL_App_QT.cpp : Definiert den Einstiegspunkt f�r die Konsolenanwendung.
//

#include <QtWidgets/QApplication>
#include "include/MainWindow.h"

int main(int argc, char* argv[])
{
  QApplication app(argc, argv);

  MainWindow glWindow;
  glWindow.resize(1280,720);

  glWindow.show();

  return app.exec();
}

