#include "ConvertSyntaxHighlighter.h"
#include <QFileInfo>
#include <QSettings>
#include <QDir>

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

  // Find things that start with a minus sign
  // QRegExp reCommand = QRegExp("\\-\\b\\[a-z]+\\b");
  QRegExp reCommand = QRegExp("(^|\\s)(-[a-z\\-]+)(\\s|$)");
  int index = text.indexOf(reCommand);
  while(index >= 0)
    {
    QString command = text.mid(reCommand.pos(2), reCommand.cap(2).length());
    if(m_CommandList.indexOf(command) >= 0)
      setFormat(reCommand.pos(2), reCommand.cap(2).length(), fmtCommand);
    index = text.indexOf(reCommand, std::max(reCommand.pos(3), index+1));
    }

  // Find things that look like filenames
  reCommand = QRegExp("(^|\\s)(\\S+)(\\s|$)");
  index = text.indexOf(reCommand);
  while(index >= 0)
    {
    QString dir = QSettings().value("working_dir", QDir::currentPath()).toString();
    QString testfile = reCommand.cap(2);
    QFileInfo qfi(dir, testfile);
    if(qfi.isFile())
      {
      setFormat(reCommand.pos(2), testfile.length(), fmtFile);
      }
    index = text.indexOf(reCommand, std::max(reCommand.pos(3), index+1));
    }
}
