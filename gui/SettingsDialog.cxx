#include "SettingsDialog.h"
#include "ui_SettingsDialog.h"
#include <QFileDialog>
#include <QStandardPaths>
#include <QSettings>

SettingsDialog::SettingsDialog(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::SettingsDialog)
{
  ui->setupUi(this);
  ui->labelInvalidViewerPath->setVisible(false);
}

SettingsDialog::~SettingsDialog()
{
  delete ui;
}

void SettingsDialog::show()
{
  QDialog::show();
  ui->inViewerPath->setText(QSettings().value("viewerPath").toString());
}

void SettingsDialog::on_btnBrowse_clicked()
{

  QString start =
      QStandardPaths::standardLocations(QStandardPaths::ApplicationsLocation).first();

  QString dir =
      QFileDialog::getOpenFileName(this, "Locate Image Viewer", start);

#ifdef __APPLE__

  // If the user selects an app bundle, we will look for an executable in the bundle
  if(dir.endsWith(".app", Qt::CaseInsensitive))
    {
    dir = QStandardPaths::findExecutable("ITK-SNAP", QStringList(dir + "/Contents/MacOS"));
    }

#endif

  ui->inViewerPath->setText(dir);
}

void SettingsDialog::on_inViewerPath_textChanged(const QString &text)
{
  // Some automatic replacements
  QString dir = text;

#ifdef __APPLE__

  // If the user selects an app bundle, we will look for an executable in the bundle
  if(dir.endsWith(".app", Qt::CaseInsensitive))
    {
    dir = QStandardPaths::findExecutable("ITK-SNAP", QStringList(dir + "/Contents/MacOS"));
    }

#endif

  // Check the path
  QFileInfo fi(dir);
  if((fi.exists() && fi.isExecutable()) || dir.length() == 0)
    ui->labelInvalidViewerPath->setVisible(false);
  else
    ui->labelInvalidViewerPath->setVisible(true);
}

void SettingsDialog::on_SettingsDialog_accepted()
{
  // When the user presses OK, we commit the path if it is valid, clear it otherwise
  if(ui->inViewerPath->text().length() && !ui->labelInvalidViewerPath->isVisible())
    {
    QSettings().setValue("viewerPath", ui->inViewerPath->text());
    }
  else
    {
    QSettings().setValue("viewerPath", QVariant());
    }
}



void SettingsDialog::on_SettingsDialog_rejected()
{
}
