#include "HistoryDialog.h"
#include "ui_HistoryDialog.h"
#include <QSettings>
#include <QDateTime>
#include <QDebug>
#include <QClipboard>

HistoryDialog::HistoryDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::HistoryDialog)
{
  ui->setupUi(this);
  m_HistoryModel = new QStandardItemModel(0, 1);
  m_HistoryModel->setHeaderData(0, Qt::Horizontal, "Command");

  m_ProxyModel = new QSortFilterProxyModel(this);
  m_ProxyModel->setSourceModel(m_HistoryModel);
  ui->tvHistory->setModel(m_ProxyModel);

  m_ProxyModel->setFilterRole(Qt::UserRole+3);
  m_ProxyModel->setSortRole(Qt::UserRole+2);
  m_ProxyModel->setDynamicSortFilter(true);
  m_ProxyModel->sort(0, Qt::DescendingOrder);

  connect(ui->tvHistory, SIGNAL(doubleClicked(QModelIndex)),
          this, SLOT(onHistoryIndexDoubleClick(QModelIndex)));

  connect(ui->inSearch, SIGNAL(textChanged(QString)),
          m_ProxyModel, SLOT(setFilterWildcard(QString)));

  // Context menu
  ui->tvHistory->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(ui->tvHistory, SIGNAL(customContextMenuRequested(QPoint)),
          this, SLOT(onHistoryContextRequest(QPoint)));
}

HistoryDialog::~HistoryDialog()
{
  delete ui;
}

void HistoryDialog::addHistoryEntry(QImage &pixmap, const QString &wd, QString &command,
                                    QVariant timestamp)
{
  if(timestamp.isNull())
    timestamp = QDateTime::currentDateTimeUtc();

  QStandardItem *item = new QStandardItem();
  // item->setData(pixmap, Qt::DecorationRole);
  item->setData(wd, Qt::UserRole+1);
  item->setData(timestamp, Qt::UserRole+2);
  item->setData(command, Qt::UserRole+3);
  item->setData(command.simplified(), Qt::DisplayRole);

  QDateTime ts = timestamp.toDateTime();

  QString tooltip = QString("Date: %2\nWorking Dir: %1").arg(wd).arg(ts.toLocalTime().toString());
  item->setData(tooltip, Qt::ToolTipRole);

  m_HistoryModel->appendRow(item);
}


void HistoryDialog::loadHistory()
{
  m_HistoryModel->clear();

  QSettings settings;
  int size = settings.beginReadArray("history");
  for(int i = 0; i < size; i++)
  {
    settings.setArrayIndex(i);
    QString command = settings.value("command").toString();
    QString wd = settings.value("workdir").toString();
    QDateTime ts = settings.value("timestamp").toDateTime();
    QImage pixmap = settings.value("pixmap").value<QImage>();

    this->addHistoryEntry(pixmap, wd, command, ts);
  }
  settings.endArray();
}

QSize HistoryDialog::sizeHint() const
{
  QSize hint = QWidget::sizeHint();
  hint.setWidth(400);
  return hint;
}

void HistoryDialog::saveHistory()
{
  QSettings settings;
  settings.beginWriteArray("history");
  for(int i = 0; i < m_HistoryModel->rowCount(); i++)
  {
    settings.setArrayIndex(i);

    QStandardItem *item = m_HistoryModel->item(i);
    settings.setValue("command", item->data(Qt::UserRole+3));
    settings.setValue("workdir", item->data(Qt::UserRole+1));
    settings.setValue("timestamp", item->data(Qt::UserRole+2));
    settings.setValue("pixmap", item->data(Qt::DecorationRole));
  }
  settings.endArray();
  settings.sync();
}


void HistoryDialog::on_btnClose_clicked()
{
  this->accept();
}

void HistoryDialog::on_btnCopy_clicked()
{
  QModelIndexList slist = ui->tvHistory->selectionModel()->selectedRows();
  if(slist.size() == 1)
    {
    int i = slist.first().row();
    QStandardItem *item = m_HistoryModel->item(i);

    QString command = item->data(Qt::UserRole+3).toString();
    emit commandCopyRequested(command);
    }
}

void HistoryDialog::onHistoryIndexDoubleClick(QModelIndex mi)
{
  int row = m_ProxyModel->mapToSource(mi).row();
  QStandardItem *item = m_HistoryModel->item(row);
  QString command = item->data(Qt::UserRole+3).toString();
  emit commandCopyRequested(command);
}

#include <QMenu>

void HistoryDialog::onHistoryContextRequest(const QPoint &point)
{
  QPoint globalPos = ui->tvHistory->viewport()->mapToGlobal(point);
  QModelIndex idx = ui->tvHistory->indexAt(point);


  if(!idx.isValid())
    return;

  QMenu myMenu;
  QAction *a_use = myMenu.addAction("Place in Editor");
  QAction *a_copy = myMenu.addAction("Copy to Clipboard");
  QAction *a_chdir = myMenu.addAction("Change to Command's Directory");
  myMenu.addSeparator();
  QAction *a_remove = myMenu.addAction("Remove from History");

  QAction* action = myMenu.exec(globalPos);
  if(!action)
    return;

  int row = m_ProxyModel->mapToSource(idx).row();
  QStandardItem *item = m_HistoryModel->item(row);

  if(action == a_chdir)
    {
    emit changeDirectoryRequested(item->data(Qt::UserRole+1).toString());
    }
  else if(action == a_remove)
    {
    m_HistoryModel->removeRow(row);
    }
  else if(action == a_use)
    {
    emit commandCopyRequested(item->data(Qt::UserRole+3).toString());
    }
  else if(action == a_copy)
    {
    QApplication::clipboard()->setText(item->data(Qt::UserRole+3).toString());
    }

}
