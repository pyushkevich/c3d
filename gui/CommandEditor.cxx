/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    CommandEditor.cxx
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

#include "CommandEditor.h"
#include <QCompleter>
#include <QAbstractItemView>
#include <QKeyEvent>
#include <QScrollBar>
#include <QDir>
#include <QUrl>
#include <QSettings>
#include <QToolTip>
#include <QDebug>
#include <QMimeData>
#include <QMouseEvent>
#include <QDateTime>
#include <QFontMetrics>
#include "Documentation.h"
#include <sstream>

CommandEditor::CommandEditor(QWidget *parent):
  QTextEdit(parent)
{
  m_fileCompleter = NULL;
  m_commandCompleter = NULL;
}

CommandEditor::~CommandEditor()
{
}

void CommandEditor::setFileCompleter(QCompleter *c)
{
  if (m_fileCompleter)
    QObject::disconnect(m_fileCompleter, 0, this, 0);

  m_fileCompleter = c;

  if (!m_fileCompleter)
    return;

  m_fileCompleter->setWidget(this);
  m_fileCompleter->setCompletionMode(QCompleter::PopupCompletion);
  QObject::connect(m_fileCompleter, SIGNAL(activated(QString)),
                   this, SLOT(insertFileCompletion(QString)));
}

QCompleter *CommandEditor::fileCompleter() const
{
  return m_fileCompleter;
}

void CommandEditor::setCommandCompleter(QCompleter *c)
{
  if (m_commandCompleter)
    QObject::disconnect(m_commandCompleter, 0, this, 0);

  m_commandCompleter = c;

  if (!m_commandCompleter)
    return;

  m_commandCompleter->setWidget(this);
  m_commandCompleter->setCompletionMode(QCompleter::PopupCompletion);
  QObject::connect(m_commandCompleter, SIGNAL(activated(QString)),
                   this, SLOT(insertCommandCompletion(QString)));
}

QCompleter *CommandEditor::commandCompleter() const
{
  return m_commandCompleter;
}

void CommandEditor::setDocumentation(Documentation *doc)
{
  m_Documentation = doc;
}


QString splitTooltip(QString text, int width) 
{ 
    QFontMetrics fm(QToolTip::font()); 
    QString result; 

    for (;;) { 
        int i = 0; 
        while (i < text.length()) { 
            if (fm.width(text.left(++i + 1)) > width) { 
                int j = text.lastIndexOf(' ', i); 
                if (j > 0) 
                    i = j; 
                result += text.left(i); 
                result += '\n'; 
                text = text.mid(i+1); 
                break; 
            } 
            else if(text[i] == '\n')
              {
              result += text.left(i+1);
              text = text.mid(i+1);
              }
        } 
        if (i >= text.length()) 
            break; 
    } 
    return result + text; 
} 

#include <iostream>
bool CommandEditor::event(QEvent *e)
{
  if (e->type() == QEvent::ToolTip)
    {
    QHelpEvent *helpEvent = static_cast<QHelpEvent *>(e);
    QTextCursor tc = this->cursorForPosition(helpEvent->pos());
    QString fn = this->filenameUnderCursor(tc);
    QString help = QString("future help for %1").arg(fn);

    // Get the c++ string
    std::string fn_std = fn.toStdString();

    // Special handling for real files
    QFileInfo fi(QDir(QSettings().value("working_dir").toString()), fn);
    if(fi.exists() && fi.isReadable())
      {
      help=QString(
        "<html><b>File: </b>%2<br>"
        "<b>Directory: </b>%1<br>"
        "<b>Last Modified:  </b>%3<br><br><br>"
        "<i>Double click to open in a 3D viewer</i></html>")
          .arg(fi.absolutePath())
          .arg(fi.fileName())
          .arg(fi.lastModified().toLocalTime().toString());
      }

    // Special handling for c3d commands
    else if(m_Documentation->GetAllCommands().find(fn_std)
      != m_Documentation->GetAllCommands().end())
      {
      std::ostringstream oss;
      std::string nodash = fn_std.substr(1);
      if(m_Documentation->PrintCommandHelp(oss, nodash.c_str()))
        {
        help = splitTooltip(QString::fromStdString(oss.str()), 600);
        }
      }

    QToolTip::showText(helpEvent->globalPos(), help);
    return true;
    }

  else
    return QTextEdit::event(e);
}

bool CommandEditor::canInsertFromMimeData(const QMimeData *source) const
{
  if(source->hasUrls())
    return true;
  else return QTextEdit::canInsertFromMimeData(source);
}

void CommandEditor::insertFromMimeData(const QMimeData *source)
{
  if(source->hasUrls())
    {
    foreach(QUrl url, source->urls())
      {
      if(url.isLocalFile())
        {
        this->textCursor().insertText(url.toLocalFile());
        this->textCursor().insertText(" ");
        }
      }
    }
  else
    {
    QTextEdit::insertFromMimeData(source);
    }
}

#include <iostream>

void CommandEditor::keyPressEvent(QKeyEvent *e)
{
  bool in_file_popup =
      m_fileCompleter && m_fileCompleter->popup() && m_fileCompleter->popup()->isVisible();
  bool in_cmd_popup =
      m_commandCompleter && m_commandCompleter->popup() && m_commandCompleter->popup()->isVisible();
  bool in_popup = in_file_popup | in_cmd_popup;

  // If a popup is open, we don't want to do anything with these keys
  if (in_popup)
    {
    // The following keys are forwarded by the completer to the widget
    switch (e->key())
      {
      case Qt::Key_Enter:
      case Qt::Key_Return:
      case Qt::Key_Escape:
      case Qt::Key_Tab:
      case Qt::Key_Backtab:
        e->ignore();
        return; // let the completer do default behavior
      default:
        break;
      }
    }

  // Command-enter - execute command
  if((e->key() == Qt::Key_Enter || e->key() == Qt::Key_Return)
     && (e->modifiers() == Qt::ControlModifier))
    {
    emit commandAccepted();
    e->ignore();
    return;
    }

  // Command-backspace - clear contents
  if((e->key() == Qt::Key_Backspace || e->key() == Qt::Key_Delete)
     && (e->modifiers() == Qt::ControlModifier))
    {
    emit clearRequested();
    e->ignore();
    return;
    }

  // Indentation. If the user pressed enter, mimic the leading spaces
  // from the previous line
  if(e->key() == Qt::Key_Enter || e->key() == Qt::Key_Return)
    {
    // Go to the beginning of the line and count the blank characters
    QTextCursor tc = textCursor();
    tc.movePosition(QTextCursor::StartOfLine, QTextCursor::MoveAnchor);

    // Indent after 'c3d ' or after spaces
    QRegExp indentme("(c[2-4]d|snap|itksnap|view)\\s*|\\s*", Qt::CaseInsensitive);
    QString leading = this->document()->find(indentme, tc).selectedText();
    QTextEdit::keyPressEvent(e);
    textCursor().insertText(QString(leading.length(),' '));
    return;
    }

  // Test for the shortcut
  bool isShortcut = (!e->modifiers() && e->key() == Qt::Key_Tab);
  if(!isShortcut && !in_popup)
    {
    QTextEdit::keyPressEvent(e);
    return;
    }

  // We also want to break out if we are not in popup mode, and the character
  // to the left of the cursor is a space.
  QTextCursor tc = textCursor();
  if(tc.position() == 0)
    {
    textCursor().insertText("    ");
    return;
    }

  tc.movePosition(QTextCursor::Left, QTextCursor::MoveAnchor);
  tc.movePosition(QTextCursor::Right, QTextCursor::KeepAnchor);
  if(tc.selectedText().indexOf(QRegExp("\\s+")) >= 0)
    {
    textCursor().insertText("    ");
    return;
    }

  // We are here because a shortcut was pressed or because we are in popup
  // mode and a key was pressed. If the latter, we want to send the key to
  // the text editor, and keep going
  if(in_popup)
    {
    QTextEdit::keyPressEvent(e);
    }

  // Check the word under the cursor
  QString userText = this->filenameUnderCursor();

  // If the text is zero length, treat it as a usual tab
  if(userText.length() == 0)
    {
    QTextEdit::keyPressEvent(e);
    return;
    }

  // What completion to do?
  if(userText.left(1) == "-")
    {
    m_commandCompleter->setCompletionPrefix(userText);
    if(m_commandCompleter->completionCount() == 1)
      {
      if(in_cmd_popup)
        m_commandCompleter->popup()->hide();
      this->insertCommandCompletion(m_commandCompleter->currentCompletion());
      }
    else
      {
      m_commandCompleter->setCompletionMode(QCompleter::PopupCompletion);
      m_commandCompleter->popup()->setCurrentIndex(
            m_commandCompleter->completionModel()->index(0, 0));

      if(!in_cmd_popup)
        {
        m_popupRect = cursorRect();
        m_popupRect.setWidth(
              m_commandCompleter->popup()->sizeHintForColumn(0)
              + m_commandCompleter->popup()->verticalScrollBar()->sizeHint().width());
        }

      m_commandCompleter->complete(m_popupRect); // popup it up!
      }
    }

  else
    {
    // Turn the completion prefix into a real filename
    QString wdir = QSettings().value("working_dir", QDir::currentPath()).toString();
    QString userPath = QDir::cleanPath(QDir(wdir).absoluteFilePath(userText));

    // If the user specifies an actual directory that exists, we want to search
    // the contents of that directory
    if(QFileInfo(userPath).isDir() &&
       (userText.right(1) == QDir::separator() || userText.right(1) == "/"))
      {
      userPath = userPath + "/";
      }

    // How many completions?
    m_fileCompleter->setCompletionPrefix(userPath);
    if(m_fileCompleter->completionCount() == 1)
      {
      m_completionFileRelative = userText;
      if(in_file_popup)
        m_fileCompleter->popup()->hide();
      this->insertFileCompletion(m_fileCompleter->currentCompletion());
      }
    else
      {
      m_completionFileRelative = userText;
      m_fileCompleter->setCompletionMode(QCompleter::PopupCompletion);
      m_fileCompleter->popup()->setCurrentIndex(m_fileCompleter->completionModel()->index(0, 0));

      if(!in_file_popup)
        {
        m_popupRect = cursorRect();
        m_popupRect.setWidth(
              m_fileCompleter->popup()->sizeHintForColumn(0)
              + m_fileCompleter->popup()->verticalScrollBar()->sizeHint().width());
        }

      m_fileCompleter->complete(m_popupRect); // popup it up!
      }
    }


  /*
  // Popup the completer
  if(userPath != m_fileCompleter->currentCompletion())
    {
    m_completionFileRelative = userText;
    m_fileCompleter->setCompletionPrefix(userPath);
    m_fileCompleter->complete();
    }
  else
    {
    m_fileCompleter->setCurrentRow(m_fileCompleter->currentRow() + 1);
    m_fileCompleter->complete();
    }
    */

  /*
  static QString eow("~!@#$%^&*()_+{}|:\"<>?,./;'[]\\-="); // end of word
  bool hasModifier = (e->modifiers() != Qt::NoModifier) && !ctrlOrShift;

  if (!isShortcut && (hasModifier || e->text().isEmpty()|| completionPrefix.length() < 3
                      || eow.contains(e->text().right(1)))) {
    m_fileCompleter->popup()->hide();
    return;
    }

  if (completionPrefix != m_fileCompleter->completionPrefix()) {
    m_fileCompleter->setCompletionPrefix(completionPrefix);
    m_fileCompleter->popup()->setCurrentIndex(m_fileCompleter->completionModel()->index(0, 0));
    }

  QRect cr = cursorRect();
  cr.setWidth(m_fileCompleter->popup()->sizeHintForColumn(0)
              + m_fileCompleter->popup()->verticalScrollBar()->sizeHint().width());
  m_fileCompleter->complete(cr); // popup it up!
  */
}

void CommandEditor::focusInEvent(QFocusEvent *e)
{
  if (m_fileCompleter)
    m_fileCompleter->setWidget(this);
  QTextEdit::focusInEvent(e);
}

void CommandEditor::mouseDoubleClickEvent(QMouseEvent *mev)
{
  QTextCursor cursor = this->cursorForPosition(mev->pos());
  QString filename = this->filenameUnderCursor(cursor);
  QFileInfo fi(QDir(QSettings().value("working_dir").toString()), filename);

  // Does the file exist?
  if(fi.exists() && fi.isFile() && fi.isReadable())
    {
    emit validFilenameClicked(fi.absoluteFilePath());
    }
  else
    QTextEdit::mousePressEvent(mev);
}

void CommandEditor::insertFileCompletion(const QString &completion)
{
  if (m_fileCompleter->widget() != this)
    return;


  // Replace the selection with the new selection
  QString prefix = m_fileCompleter->completionPrefix();
  int extra = completion.length() - prefix.length();
  QString newtext = m_completionFileRelative + completion.right(extra);
  m_activeCompletion.removeSelectedText();
  m_activeCompletion.insertText(newtext);
}

void CommandEditor::insertCommandCompletion(const QString &completion)
{
  if (m_commandCompleter->widget() != this)
    return;


  // Replace the selection with the new selection
  QString prefix = m_commandCompleter->completionPrefix();
  int extra = completion.length() - prefix.length();
  QString newtext = prefix + completion.right(extra);
  m_activeCompletion.removeSelectedText();
  m_activeCompletion.insertText(newtext);
}

QString CommandEditor::filenameUnderCursor(QTextCursor tc)
{
  // Find the text around the cursor that constitutes a filename

  // Get the cursor
  if(tc.isNull())
    tc = textCursor();

  // Find the beginning and end of the current line
  QTextCursor cline = tc;
  cline.movePosition(QTextCursor::StartOfLine, QTextCursor::MoveAnchor);
  cline.movePosition(QTextCursor::EndOfLine, QTextCursor::KeepAnchor);

  // Search backwards to find the first whitespace character
  QTextCursor prev =
      this->document()->find(QRegExp("\\s"), tc, QTextDocument::FindBackward);

  // Search forward for next whitespace character
  QTextCursor next =
      this->document()->find(QRegExp("\\s"), tc);

  // Set the position
  int a = (prev.isNull() || prev.anchor() < cline.anchor())
      ? cline.anchor() : prev.anchor() + 1;

  int p = (next.isNull() || next.position() > cline.position())
      ? cline.position() : next.position() - 1;

  tc.setPosition(a, QTextCursor::MoveAnchor);
  tc.setPosition(p, QTextCursor::KeepAnchor);

  // Save the active completion
  m_activeCompletion = tc;

  return tc.selectedText();
}


