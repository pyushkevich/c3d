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
