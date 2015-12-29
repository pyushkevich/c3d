/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    MainWindow.h
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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QProcess>
#include <QFileSystemModel>

namespace Ui {
class MainWindow;
}

class QListView;
class HistoryDialog;
class SettingsDialog;
class Documentation;

class MainWindow : public QMainWindow
{
  Q_OBJECT
  
public:
  explicit MainWindow(QWidget *parent = 0);
  ~MainWindow();

signals:

  void highlightingInvalidated();

public slots:

  void onCommandReceive(QString command);
  void onImageViewRequested(QString filename);
  void onCommandActionTriggered();

private slots:
  void on_btnChangeDir_clicked();
  void onWorkingDirectoryChanged(const QString &dir);

  void on_btnRun_clicked();

  void onProcessFinished(int rc);
  void onProcessStandardOutput();
  void onProcessStandardError();
  void onProcessFailed(QProcess::ProcessError errorCode);

  void onFileListDoubleClick(const QModelIndex &index);

  void on_btnClear_clicked();

  void on_actionPreferences_triggered();

  void on_actionC3D_Manual_triggered();

protected:

private:
  Ui::MainWindow *ui;
  QFileSystemModel *m_FileSystemModel;
  QListView *m_CurrentDirView;
  void setupSearchPaths();

  HistoryDialog *m_History;
  SettingsDialog *m_Settings;
  Documentation *m_Documentation;
};

#endif // MAINWINDOW_H
