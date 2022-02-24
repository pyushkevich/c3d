/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ConvertSyntaxHighlighter.cxx
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

#include "ConvertSyntaxHighlighter.h"
#include <QFileInfo>
#include <QSettings>
#include <QDir>
#include <QRegularExpression>

ConvertSyntaxHighlighter::ConvertSyntaxHighlighter(QTextDocument *parent) :
  QSyntaxHighlighter(parent)
{
}

void ConvertSyntaxHighlighter::setCommandList(const QStringList &cl)
{
  m_CommandList = cl;
}

void ConvertSyntaxHighlighter::highlightBlock(const QString &text)
{
  // Format for commands
  QTextCharFormat fmtCommand;
  fmtCommand.setFontWeight(QFont::Bold);
  fmtCommand.setForeground(QColor("green"));

  // Format for filenames
  QTextCharFormat fmtFile;
  fmtFile.setUnderlineStyle(QTextCharFormat::SingleUnderline);
  fmtFile.setForeground(QColor("blue"));

  // Format for filenames
  QTextCharFormat fmtExe;
  fmtExe.setFontWeight(QFont::Bold);
  fmtExe.setForeground(QColor("darkred"));

  // Find things that start with a minus sign
  // QRegExp reCommand = QRegExp("\\-\\b\\[a-z]+\\b");
  QRegularExpressionMatch match;
  QRegularExpression reCommand1("(^|\\s)(-[a-z\\-]+)(\\s|$)");
  int index = text.indexOf(reCommand1, 0, &match);
  while(index >= 0)
    {
    QString command = text.mid(match.capturedStart(2), match.capturedLength(2));
    if(m_CommandList.indexOf(command) >= 0)
      setFormat(match.capturedStart(2), match.capturedLength(2), fmtCommand);
    index = text.indexOf(reCommand1, std::max((int) match.capturedStart(3), index+1), &match);
    }

  // Find things that look like filenames
  QRegularExpression reCommand2("(^|\\s)(\\S+)(\\s|$)");
  index = text.indexOf(reCommand2, 0, &match);
  while(index >= 0)
    {
    QString dir = QSettings().value("working_dir", QDir::currentPath()).toString();
    QString testfile = match.captured(2);
    QFileInfo qfi(dir, testfile);
    if(qfi.isFile())
      {
      setFormat(match.capturedStart(2), testfile.length(), fmtFile);
      }
    index = text.indexOf(reCommand2, std::max((int) match.capturedStart(3), index+1), &match);
    }

  index = text.indexOf(QRegularExpression("(^|\\s)(c[2-4]d|snap|itksnap|view)(\\s|$)"), 0, &match);
  if(index >= 0)
    {
    setFormat(match.capturedStart(2), match.capturedLength(2), fmtExe);
    }
}
