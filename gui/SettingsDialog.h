#ifndef SETTINGSDIALOG_H
#define SETTINGSDIALOG_H

#include <QDialog>

namespace Ui {
  class SettingsDialog;
}

class SettingsDialog : public QDialog
{
  Q_OBJECT

public:
  explicit SettingsDialog(QWidget *parent = 0);
  ~SettingsDialog();

public slots:

  virtual void show();

private slots:
  void on_btnBrowse_clicked();

  void on_inViewerPath_textChanged(const QString &arg1);

  void on_SettingsDialog_accepted();

  void on_SettingsDialog_rejected();

private:
  Ui::SettingsDialog *ui;
};

#endif // SETTINGSDIALOG_H
