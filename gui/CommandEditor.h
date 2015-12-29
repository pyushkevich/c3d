/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    CommandEditor.h
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

#ifndef COMMANDEDITOR_H
#define COMMANDEDITOR_H

#include <QTextEdit>
#include <QTextCursor>

class QCompleter;
class Documentation;

class CommandEditor : public QTextEdit
{
  Q_OBJECT
public:
  explicit CommandEditor(QWidget *parent = 0);
  
  ~CommandEditor();

  void setFileCompleter(QCompleter *c);
  QCompleter *fileCompleter() const;

  void setCommandCompleter(QCompleter *c);
  QCompleter *commandCompleter() const;

  void setDocumentation(Documentation *doc);

  bool event(QEvent *e);

  bool canInsertFromMimeData(const QMimeData *source) const;
  void insertFromMimeData(const QMimeData *source);

signals:

  void commandAccepted();
  void clearRequested();
  void validFilenameClicked(QString);

protected:
  void keyPressEvent(QKeyEvent *e);
  void focusInEvent(QFocusEvent *e);
  void mouseDoubleClickEvent(QMouseEvent *mev);

private slots:
  void insertFileCompletion(const QString &completion);
  void insertCommandCompletion(const QString &completion);

private:
  QString filenameUnderCursor(QTextCursor tc = QTextCursor());

private:
  QCompleter *m_fileCompleter, *m_commandCompleter;
  QTextCursor m_activeCompletion;
  QString m_completionFileRelative;
  QRect m_popupRect;

  QStringList m_CommandList;

  /** Documentation system */
  Documentation *m_Documentation;

signals:
  
public slots:
  
};

#endif // COMMANDEDITOR_H
