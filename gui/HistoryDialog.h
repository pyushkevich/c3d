#ifndef HISTORYDIALOG_H
#define HISTORYDIALOG_H

#include <QDialog>
#include <QStandardItemModel>

namespace Ui {
class HistoryDialog;
}

class HistoryDialog : public QDialog
{
  Q_OBJECT

signals:

  void commandCopyRequested(QString);

public:
  explicit HistoryDialog(QWidget *parent = 0);
  ~HistoryDialog();

  void addHistoryEntry(QImage &pixmap, const QString &wd, QString &command,
                       QVariant timestamp = QVariant());

  void saveHistory();
  void loadHistory();
protected:

  QStandardItemModel *m_HistoryModel;

private slots:

  void on_btnClose_clicked();

  void on_btnCopy_clicked();

  void onHistoryIndexDoubleClick(QModelIndex mi);

private:
  Ui::HistoryDialog *ui;
};

#endif // HISTORYDIALOG_H
