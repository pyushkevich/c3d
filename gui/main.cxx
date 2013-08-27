#include <QApplication>
#include <QSettings>
#include "MainWindow.h"

int main(int argc, char *argv[])
{
  QCoreApplication::setOrganizationName("org.itksnap");
  QCoreApplication::setOrganizationDomain("itksnap.org");
  QCoreApplication::setApplicationName("Convert3DGUI");

  QApplication app(argc, argv);
  MainWindow *mwin = new MainWindow();
  mwin->show();
  app.exec();
}
