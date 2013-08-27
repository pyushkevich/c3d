#include "MainWindow.h"
#include "ui_MainWindow.h"

#include "ConvertSyntaxHighlighter.h"
#include <QDir>
#include <QFileSystemModel>
#include <QStringListModel>
#include <QCompleter>
#include <QFontMetrics>
#include <QSettings>
#include <QFileDialog>
#include <QProcess>
#include <QHelpEvent>
#include "ConvertImageND.h"

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow)
{
  ui->setupUi(this);

  // Get the list of available commands
  std::ostringstream oss;
  ImageConverter<double,3>::PrintCommandListing(oss);
  QString cmdliststr(oss.str().c_str());
  QStringList cmdlist = cmdliststr.split(QRegExp("(,|\\s)"), QString::SkipEmptyParts);

  // Set the font in the editors
  QFont font;
  font.setFamily("Courier");
  font.setFixedPitch(true);
  font.setPointSize(12);

  ui->teCommand->setFont(font);
  ui->teOutput->setFont(font);

  QFontMetrics metrics(font);
  ui->teCommand->setTabStopWidth(4 * metrics.width(' '));

  // Setup the highlighter
  ConvertSyntaxHighlighter *shl =
    new ConvertSyntaxHighlighter(ui->teCommand->document());
  shl->setCommandList(cmdlist);

  // Setup the filename completer in the editor
  QCompleter *fileCompleter = new QCompleter(this);
  m_FileSystemModel = new QFileSystemModel(fileCompleter);
  fileCompleter->setModel(m_FileSystemModel);
  ui->teCommand->setFileCompleter(fileCompleter);

  // Set the command completer
  QCompleter *commandCompleter = new QCompleter(this);
  QStringListModel *cmdlistmodel = new QStringListModel(commandCompleter);
  cmdlistmodel->setStringList(cmdlist);
  commandCompleter->setModel(cmdlistmodel);
  ui->teCommand->setCommandCompleter(commandCompleter);

  // Set the directory info
  QSettings settings;
  QString dir = settings.value("working_dir", QDir::currentPath()).toString();
  this->onWorkingDirectoryChanged(dir);
}

MainWindow::~MainWindow()
{
  delete ui;
}


void MainWindow::on_btnChangeDir_clicked()
{
  QString dir =
      QFileDialog::getExistingDirectory(this, tr("Choose Working Directory"),
                                        ui->inWorkDir->text(), 0);
  if(QFileInfo(dir).isDir())
    {
    QSettings().setValue("working_dir", dir);
    onWorkingDirectoryChanged(dir);
    }
}

void MainWindow::onWorkingDirectoryChanged(const QString &dir)
{
  ui->inWorkDir->setText(dir);

  m_FileSystemModel->setRootPath(dir);

  QModelIndex index = m_FileSystemModel->index(dir);
  ui->teCommand->fileCompleter()->setCurrentRow(index.row());
}

void MainWindow::on_btnRun_clicked()
{
  QString program = "c3d";
  QString argstr = ui->teCommand->document()->toPlainText();
  QStringList args = argstr.split(QRegExp("(\\s|^|$)+"), QString::SkipEmptyParts);

  // Create the C3D process
  QProcess *myProcess = new QProcess(this);
  myProcess->setWorkingDirectory(ui->inWorkDir->text());
  // myProcess->setProcessChannelMode(QProcess::MergedChannels);

  // Take the output from the process
  connect(myProcess, SIGNAL(finished(int)), this, SLOT(onProcessFinished(int)));
  connect(myProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(onProcessStandardOutput()));
  connect(myProcess, SIGNAL(readyReadStandardError()), this, SLOT(onProcessStandardError()));

  // Clear the output
  ui->teOutput->clear();

  // Disable the execute button
  ui->btnRun->setEnabled(false);

  // Run the process
  myProcess->start(program, args);
}

void MainWindow::onProcessFinished(int rc)
{
  // Disable the execute button
  ui->btnRun->setEnabled(true);
}

void MainWindow::onProcessStandardOutput()
{
  QProcess *proc = qobject_cast<QProcess *>(sender());
  if(proc)
    {
    QString text(proc->readAllStandardOutput());
    QTextCursor tc = ui->teOutput->textCursor();
    QTextCharFormat fmt;
    fmt.setForeground(QColor("black"));
    tc.setCharFormat(fmt);
    tc.insertText(text);
    }
}

void MainWindow::onProcessStandardError()
{
  QProcess *proc = qobject_cast<QProcess *>(sender());
  if(proc)
    {
    QString text(proc->readAllStandardError());
    QTextCursor tc = ui->teOutput->textCursor();
    QTextCharFormat fmt;
    fmt.setForeground(QColor("red"));
    tc.setCharFormat(fmt);
    tc.insertText(text);
    }
}
