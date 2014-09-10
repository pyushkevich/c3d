#ifndef HISTORYDIALOG_H
#define HISTORYDIALOG_H

#include <QDialog>
#include <QStandardItemModel>
#include <QSortFilterProxyModel>

namespace Ui {
class HistoryDialog;
}

class HistoryDialog : public QDialog
{
  Q_OBJECT

signals:

  void commandCopyRequested(QString);
  void changeDirectoryRequested(QString);

public:
  explicit HistoryDialog(QWidget *parent = 0);
  ~HistoryDialog();

  void addHistoryEntry(QImage &pixmap, const QString &wd, QString &command,
                       QVariant timestamp = QVariant());

  void saveHistory();
  void loadHistory();

  virtual QSize sizeHint() const;
protected:

  QStandardItemModel *m_HistoryModel;
  QSortFilterProxyModel *m_ProxyModel;

private slots:

  void on_btnClose_clicked();

  void on_btnCopy_clicked();

  void onHistoryIndexDoubleClick(QModelIndex mi);

  void onHistoryContextRequest(const QPoint &point);

private:
  Ui::HistoryDialog *ui;
};

#endif // HISTORYDIALOG_H
