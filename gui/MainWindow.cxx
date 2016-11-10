/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    MainWindow.cxx
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
#include <QDockWidget>
#include <QListView>
#include <QDebug>
#include <QMessageBox>
#include <QShortcut>
#include <HistoryDialog.h>
#include <SettingsDialog.h>
#include "ConvertImageND.h"

#include "Documentation.h"

// Markdown documentation string generated at compile-time
// this looks a little odd, but works - the include file contains raw bytes
unsigned char c3d_md[] = {
  #include "markdown_docs.h"
  0x00 
};

class C3DFileSystemModel : public QFileSystemModel
{
public:
  C3DFileSystemModel(QObject *parent) : QFileSystemModel(parent) {}

  virtual Qt::DropActions supportedDragActions() const 
    { return Qt::CopyAction; }
};

void MainWindow::setupSearchPaths()
{
  // Search for the c3d executables
  QString appPath = QCoreApplication::applicationDirPath();

#ifdef __APPLE__
  QDir::addSearchPath("c3d", QDir(appPath + "/../bin").absolutePath());
  QDir::addSearchPath("c3d", appPath);
  QDir::addSearchPath("c3d", QDir::currentPath());
#else
  QDir::addSearchPath("c3d", appPath);
#endif
}

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow)
{
  ui->setupUi(this);

  // Set up the search paths
  this->setupSearchPaths();

  // Parse the markdown documentation
  m_Documentation = new Documentation(c3d_md);

  // Get the list of available commands
  const std::set<std::string> &cmdset = m_Documentation->GetAllCommands();
  QStringList cmdlist;
  for(std::set<std::string>::const_iterator it = cmdset.begin(); it != cmdset.end(); ++it)
    cmdlist.push_back(QString::fromStdString(*it));

  // Pass the documentation to the editor
  ui->teCommand->setDocumentation(m_Documentation);

  // Set the font in the editors
  QFont font;
  font.setFamily("Courier");
  font.setFixedPitch(true);
  font.setPixelSize(12);

  ui->teCommand->setFont(font);

  QFontMetrics metrics(font);
  ui->teCommand->setTabStopWidth(4 * metrics.width(' '));

  // Viewer
  connect(ui->teCommand, SIGNAL(validFilenameClicked(QString)), this, SLOT(onImageViewRequested(QString)));

  // Enable drag and drop to command editor
  ui->teCommand->viewport()->setAcceptDrops(true);

  // Disable wrapping in the output
  ui->teOutput->setFont(font);
  ui->teOutput->setWordWrapMode(QTextOption::NoWrap);

  // Setup the highlighter
  ConvertSyntaxHighlighter *shl =
    new ConvertSyntaxHighlighter(ui->teCommand->document());
  shl->setCommandList(cmdlist);

  // Highlighted needs to be aware of invalidations
  connect(this, SIGNAL(highlightingInvalidated()),
          shl, SLOT(rehighlight()));

  // Setup the filename completer in the editor
  QCompleter *fileCompleter = new QCompleter(this);
  m_FileSystemModel = new C3DFileSystemModel(fileCompleter);
  m_FileSystemModel->setFilter(QDir::AllDirs | QDir::Files | QDir::NoDot);
  fileCompleter->setModel(m_FileSystemModel);
  ui->teCommand->setFileCompleter(fileCompleter);

  // Set the command completer
  QCompleter *commandCompleter = new QCompleter(this);
  QStringListModel *cmdlistmodel = new QStringListModel(commandCompleter);
  cmdlistmodel->setStringList(cmdlist);
  commandCompleter->setModel(cmdlistmodel);
  ui->teCommand->setCommandCompleter(commandCompleter);

  // Add a dock for the file list
  QDockWidget *dock = new QDockWidget(this);
  dock->setAllowedAreas(Qt::RightDockWidgetArea | Qt::LeftDockWidgetArea);
  addDockWidget(Qt::RightDockWidgetArea, dock);

  // Create a list widget with a file model
  m_CurrentDirView = new QListView(this);
  m_CurrentDirView->setModel(m_FileSystemModel);
  m_CurrentDirView->setDragEnabled(true);
  dock->setWidget(m_CurrentDirView);
  dock->setWindowTitle("Working Directory");

  connect(m_CurrentDirView, SIGNAL(doubleClicked(QModelIndex)),
          this, SLOT(onFileListDoubleClick(QModelIndex)));

  // Set the directory info
  QSettings settings;
  QString dir = settings.value("working_dir", QDir::currentPath()).toString();
  this->onWorkingDirectoryChanged(dir);

  // Set up a shortcut for execute command (what a pain!)
  connect(ui->teCommand, SIGNAL(commandAccepted()), ui->btnRun, SLOT(click()));
  connect(ui->teCommand, SIGNAL(clearRequested()), ui->btnClear, SLOT(click()));

  // Set up the history dialog
  m_History = new HistoryDialog(this);
  connect(m_History, SIGNAL(commandCopyRequested(QString)),
          this, SLOT(onCommandReceive(QString)));
  connect(m_History, SIGNAL(changeDirectoryRequested(QString)),
          this, SLOT(onWorkingDirectoryChanged(QString)));

  // Load the history from settings
  m_History->loadHistory();

  // Add another dock for the history
  QDockWidget *dock2 = new QDockWidget(this);
  dock2->setAllowedAreas(Qt::RightDockWidgetArea | Qt::LeftDockWidgetArea);
  addDockWidget(Qt::RightDockWidgetArea, dock2);
  dock2->setWidget(m_History);
  dock2->setWindowTitle("Command History");

  // Settings
  m_Settings = new SettingsDialog(this);

  // Build the commands list
  const std::vector<Documentation::Category> &cat = m_Documentation->GetCategories();
  for(int i = 0; i < cat.size(); i++)
    {
    // Create a submenu for these commands
    QString title = QString("Category: ") + QString::fromStdString(cat[i].Title);
    QMenu *mCat = new QMenu(title, this);

    // Add all commands to this submenu
    QList<QAction *> actions;

    for(int j = 0; j < cat[i].Commands.size(); j++)
      {
      const Documentation::CommandDoc &cd = cat[i].Commands[j];
      QString mainAlias = QString::fromStdString(cd.Aliases.front());

      QAction *action = new QAction(this);
      action->setData(mainAlias);
      action->setText(QString::fromStdString(cd.Title));
      action->setToolTip(QString::fromStdString(cd.ShortDesc));
      connect(action, SIGNAL(triggered()), this, SLOT(onCommandActionTriggered()));
      actions.push_back(action);
      }

    mCat->addActions(actions);
    ui->menuCommands->addMenu(mCat);
    }
}

#include <QFileSystemModel>

MainWindow::~MainWindow()
{
  m_History->saveHistory();
  delete ui;
  delete m_Documentation;
}

void MainWindow::onCommandReceive(QString command)
{
  ui->teCommand->setPlainText(command);
  ui->teOutput->clear();
}


void MainWindow::on_btnChangeDir_clicked()
{
  QString dir =
      QFileDialog::getExistingDirectory(this, tr("Choose Working Directory"),
                                        ui->inWorkDir->text(), 0);
  if(QFileInfo(dir).isDir())
    {
    onWorkingDirectoryChanged(dir);
    }
}

void MainWindow::onWorkingDirectoryChanged(const QString &dir)
{
  QSettings().setValue("working_dir", dir);

  ui->inWorkDir->setText(dir);

  m_FileSystemModel->setRootPath(dir);

  QModelIndex index = m_FileSystemModel->index(dir);
  ui->teCommand->fileCompleter()->setCurrentRow(index.row());

  m_CurrentDirView->setRootIndex(index);

  emit highlightingInvalidated();


}

void MainWindow::onImageViewRequested(QString filename)
{
  // Get the viewer path
  QString viewerPath = QSettings().value("viewerPath").toString();
  QFileInfo fiViewer(viewerPath);

  if(!fiViewer.exists() || !fiViewer.isExecutable())
    {
    QMessageBox::warning(this, "Cannot Launch Viewer", "Can not launch viewer. Use the Preferences dialog to set viewer program");
    return;
    }

  // Create the C3D process
  QProcess *myProcess = new QProcess(this);
  myProcess->setWorkingDirectory(ui->inWorkDir->text());
  QStringList args;
  args.push_back(filename);
  myProcess->start(fiViewer.absoluteFilePath(), args);
}

void MainWindow::onCommandActionTriggered()
{
  QAction *action = (QAction *) this->sender();
  QString command = action->data().toString();

  QTextCursor tc = ui->teCommand->textCursor();
  tc.insertText(command);
  if(!ui->teCommand->document()->characterAt(tc.position()+1).isSpace())
    tc.insertText(" ");
  ui->teCommand->setFocus();
}

void MainWindow::on_btnRun_clicked()
{
  // Contents of the editor
  QString contents = ui->teCommand->document()->toPlainText();

  // Split contents into words
  QStringList words = contents.split(QRegExp("(\\s|^|$)+"), QString::SkipEmptyParts);
  if(words.size() == 0)
    return;

  // Devise the command to call
  QString command = "c3d";
  if(words[0].toLower() == "c3d" || words[0].toLower() == "c2d" || words[0].toLower() == "c4d")
    {
    command = words[0];
    words.removeFirst();
    }
  else if(words[0].toLower() == "snap" || words[0].toLower() == "itksnap"|| words[0].toLower() == "view")
    {
    command = "snap";
    words.removeFirst();
    }

  // Find the command
  QFileInfo fiCommand;
  if(command == "snap")
    {
    fiCommand.setFile(QSettings().value("viewerPath").toString());
    }
  else
    {
#ifdef _WIN32
    QFile fCommand(QString("c3d:%1.exe").arg(command));
#else
    QFile fCommand(QString("c3d:%1").arg(command));
#endif
    fiCommand.setFile(fCommand);
    }

  // Verify that C3D command exists
  if(!fiCommand.exists())
    {
    QMessageBox::warning(this, "Convert3DGUI",
                         QString("Unable to find executable %1").arg(fiCommand.fileName()));
    return;
    }

  // Create the C3D process
  QProcess *myProcess = new QProcess(this);
  myProcess->setWorkingDirectory(ui->inWorkDir->text());
  // myProcess->setProcessChannelMode(QProcess::MergedChannels);

  // Take the output from the process
  connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this , SLOT(onProcessFailed(QProcess::ProcessError)));
  connect(myProcess, SIGNAL(finished(int)), this, SLOT(onProcessFinished(int)));
  connect(myProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(onProcessStandardOutput()));
  connect(myProcess, SIGNAL(readyReadStandardError()), this, SLOT(onProcessStandardError()));

  // Clear the output
  ui->teOutput->clear();

  // Disable the execute button
  ui->btnRun->setEnabled(false);

  // Run the process
  myProcess->start(fiCommand.absoluteFilePath(), words);
}

#include <QImage>
#include <QPainter>
#include <QRgb>

void MainWindow::onProcessFinished(int rc)
{
  // Save the command to the command history
  QString wd = QSettings().value("working_dir").toString();
  QString command = ui->teCommand->document()->toPlainText();

  // Focus off!

  // Grab a picture of the screen
  QSize docsize(ui->teCommand->document()->idealWidth() + 2, ui->teCommand->document()->size().height() + 2);
  QImage image(docsize, QImage::Format_ARGB32);

  if(rc == 0)
    image.fill(Qt::white);
  else
    image.fill(QColor(0xff,0xee,0xee,0xff));

  QPaintDevice *old = ui->teCommand->document()->documentLayout()->paintDevice();
  ui->teCommand->document()->documentLayout()->setPaintDevice(&image);

  QPainter p(&image);
  ui->teCommand->document()->drawContents(&p);
  p.end();

  ui->teCommand->document()->documentLayout()->setPaintDevice(old);

  m_History->addHistoryEntry(image, wd, command);

  // Disable the execute button
  ui->btnRun->setEnabled(true);

  // Update highlighting
  emit highlightingInvalidated();

  // A message for the user
  QTextCursor tc = ui->teOutput->textCursor();
  QTextCharFormat fmt;
  QProcess *process = qobject_cast<QProcess *>(this->sender());

  if(rc != 0)
    {
    fmt.setForeground(QColor("darkred"));
    tc.setCharFormat(fmt);
    tc.insertText(QString("Command %1 returned with error code %2").arg(process->program()).arg(rc));
    }



}

void MainWindow::onProcessFailed(QProcess::ProcessError errorCode)
{
  QMessageBox::warning(this, "Convert3DGUI",
                       QString("Convert3D process failed with error %1").arg(errorCode));

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

void MainWindow::onFileListDoubleClick(const QModelIndex &index)
{
  if(m_FileSystemModel->isDir(index))
    {
    this->onWorkingDirectoryChanged(m_FileSystemModel->fileInfo(index).absoluteFilePath());
    return;
    }
  else
    {
    QString filename = m_FileSystemModel->data(index, Qt::DisplayRole).toString();
    QTextCursor tc = ui->teCommand->textCursor();
    tc.insertText(filename);
    if(!ui->teCommand->document()->characterAt(tc.position()+1).isSpace())
      tc.insertText(" ");
    m_CurrentDirView->clearFocus();
    ui->teCommand->setFocus();
  }
}

void MainWindow::on_btnClear_clicked()
{
  // Clear the document
  ui->teCommand->document()->clear();
  ui->teOutput->document()->clear();
}

void MainWindow::on_actionPreferences_triggered()
{
  m_Settings->show();
  m_Settings->raise();
}

#include <QDesktopServices>

void MainWindow::on_actionC3D_Manual_triggered()
{
  QDesktopServices::openUrl(QUrl("http://www.itksnap.org/pmwiki/pmwiki.php?n=Convert3D.Documentation"));
}
