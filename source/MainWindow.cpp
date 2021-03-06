#include "include/MainWindow.h"

#include <fstream>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QFileDialog>

//#include "moc_MainWindow.hxx"

MainWindow::MainWindow()
{
  m_glWidget = new GLwidget();
  setCentralWidget(m_glWidget);

  m_fileMenu = menuBar()->addMenu(tr("&File"));
  m_fileMenu->addAction("open",this,SLOT(openFile()));

  m_viewMenu = menuBar()->addMenu(tr("&View"));
  m_viewMenu->addAction("projection", this, SLOT(changeProjection()));
}

void MainWindow::openFile()
{
  const QString extfilter = ("Pointclouds (*.xyz *.xyzc)");
  QStringList filenames = QFileDialog::getOpenFileNames(this, "Open File", "data", extfilter, 0);

  if (filenames.empty()) return;

  std::vector<Point3d> points;
  int startDim = 0; // start Dimension for KDTree

  clock_t begin = clock();
  loadFileXYZ(filenames.front().toLocal8Bit(), points);
  clock_t end = clock();
  std::cout << "Time needed to load data: " << double(end - begin) / CLOCKS_PER_SEC << "s" << std::endl;


  begin = clock();
  KDTree dataTree = KDTree(points, startDim);
  end = clock();
  std::cout << "Time needed to build kdTree: " << double(end - begin) / CLOCKS_PER_SEC << "s" << std::endl;

  //KDTree - Abfrage
  //----------------------------------------------------------------------------
  Point3d abfragePoint = points[points.size() / 2];
  Point3d S = m_glWidget->maxPoint() - m_glWidget->minPoint();
  //TODO: Why is S = 0 0 0 ?
  double abfrageLaenge = S.x * 0.25;
  //TODO: WHY????
  //double abfrageLaenge = 1;

  std::cout << "abgfrageLaenge: " << abfrageLaenge << std::endl
            << "S:" << S.x << " " << S.y << " " << S.z << std::endl;

  begin = clock();
  std::vector<Point3d> res = dataTree.abfrage(abfrageLaenge, abfragePoint, startDim);
  end = clock();
  std::cout << "Time time needed to process subsection: " << double(end - begin) / CLOCKS_PER_SEC << "s" << std::endl;
  //----------------------------------------------------------------------------

  m_glWidget->setPoints(points);
  m_glWidget->drawPoints(res, 4, 255, 0, 0);
}

void MainWindow::changeProjection()
{
    printf("change projection\n");

    m_glWidget->makeCurrent();

    if (m_glWidget->camera().usesPerspectiveProjection())
    m_glWidget->camera().usePerspectiveProjection(false);
    else
    m_glWidget->camera().usePerspectiveProjection(true);

    m_glWidget->update();
}

//Here is the implementation of our file reader
void MainWindow::loadFileXYZ(const char* filename, std::vector<Point3d>& points)
{
    points.clear();
    std::ifstream file(filename);
    if(!file)
    {
        std::cout << "File not found\t " << filename << " is unavailable." << std::endl;
        return;
    }

    std::cout << "reading file: " << filename << std::endl;
    double a,b,c;
    while(file >> a >> b >> c)
        points.push_back(Point3d(a,b,c));
//        points.emplace_back(a,b,c);// equal to the cmd above

    size_t numberOfPoints = points.size();
    std::cout << "reading finished: " << numberOfPoints << " points have be read" << std::endl;
}
