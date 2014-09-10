/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    main.cxx
  Language:  C++
  Website:   itksnap.org/c3d
  Copyright (c) 2014 Paul A. Yushkevich
  
  This file is part of C3D, a command-line companion tool to ITK-SNAP

  C3D is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

#include <QApplication>
#include <QSettings>
#include <QStandardPaths>
#include "MainWindow.h"

void findViewer()
{
  // Name list to look through
  QStringList appNameList;

  // Start with the one saved in the system
  appNameList.append(QSettings().value("viewerPath").toString());

#ifdef __APPLE__

  // This code tries to find ITK-SNAP
  QStringList appDirList = QStandardPaths::standardLocations(QStandardPaths::ApplicationsLocation);
  foreach (QString path, appDirList)
    {
    appNameList.append(QString("%1/ITK-SNAP.app/Contents/MacOS/ITK-SNAP").arg(path));
    appNameList.append(QString("%1/ITK-SNAP.app/Contents/MacOS/InsightSNAP").arg(path));
    }

  appNameList.append("/usr/local/bin/itksnap");
  appNameList.append("/usr/local/bin/snap");
  appNameList.append("/usr/bin/itksnap");
  appNameList.append("/usr/bin/snap");

#elif _WIN32

  // Look in Program Files
  QStringList appDirList;
  appDirList.append(qgetenv("PROGRAMFILES").constData());
  appDirList.append(qgetenv("PROGRAMFILES(x86)").constData());

  foreach (QString path, appDirList)
    {
    QDir dir(path);
    if(dir.exists())
      {
      dir.setNameFilters(QStringList("ITK-SNAP*"));
      dir.setSorting(QDir::Time | QDir::Reversed);
      QStringList matches = dir.entryList();
      foreach(QString prog, matches)
        {
        appNameList.append(QString("%1\\%2\\bin\\ITK-SNAP.exe").arg(path).arg(prog));
        }
      }
    }

#else

  // This 

#endif

  // Search through all options
  foreach (QString prog, appNameList)
    {
    QFileInfo fi(prog);
    if(fi.exists() && fi.isExecutable())
      {
      QSettings().setValue("viewerPath", prog);
      break;
      }
    }
}

int main(int argc, char *argv[])
{
  QCoreApplication::setOrganizationName("org.itksnap");
  QCoreApplication::setOrganizationDomain("itksnap.org");
  QCoreApplication::setApplicationName("Convert3DGUI");

  // Find ITK-SNAP
  findViewer();

  QApplication app(argc, argv);
  MainWindow *mwin = new MainWindow();
  mwin->show();
  app.exec();

  delete mwin;
}
