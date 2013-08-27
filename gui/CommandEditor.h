#ifndef COMMANDEDITOR_H
#define COMMANDEDITOR_H

#include <QTextEdit>
#include <QTextCursor>

class QCompleter;

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

  void setCommandList(const QStringList &cl);

  bool event(QEvent *e);

protected:
  void keyPressEvent(QKeyEvent *e);
  void focusInEvent(QFocusEvent *e);

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

signals:
  
public slots:
  
};

#endif // COMMANDEDITOR_H
