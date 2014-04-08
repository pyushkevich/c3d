/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ConvertSyntaxHighlighter.h
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

#ifndef CONVERTSYNTAXHIGHLIGHTER_H
#define CONVERTSYNTAXHIGHLIGHTER_H

#include <QSyntaxHighlighter>

class ConvertSyntaxHighlighter : public QSyntaxHighlighter
{
  Q_OBJECT
public:
  explicit ConvertSyntaxHighlighter(QTextDocument *parent = 0);

  void setCommandList(const QStringList &cl);

  virtual void highlightBlock(const QString &text);
  
signals:
  
public slots:

private:
  QStringList m_CommandList;
  
};

#endif // CONVERTSYNTAXHIGHLIGHTER_H
