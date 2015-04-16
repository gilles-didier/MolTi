/*
    'MolTi' and 'molti-console' detects communities from multiplex networks / 'bonf' computes q-values of annotations enrichment of communities / 'test' simulates random multiplex to test community detection approaches
    Copyright (C) 2015  Gilles DIDIER

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#include "moltiwindow.h"
#include "ui_moltiwindow.h"
#include "qitemclassdelegate.h"
#include "qitemdoubledelegate.h"
#include "qitemintdelegate.h"
#include "qitembigstringdelegate.h"
#include <QtWidgets>

#define HELP_PAGE "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">\n<html>\n<head>\n	<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\"/>\n	<title></title>\n	<meta name=\"generator\" content=\"LibreOffice 4.3.3.2 (Linux)\"/>\n	<meta name=\"created\" content=\"2015-03-31T00:00:00\"/>\n	<meta name=\"changed\" content=\"2015-04-01T11:20:48.923424128\"/>\n	<style type=\"text/css\">\n		@page { margin: 2cm }\n		p { margin-bottom: 0.25cm; line-height: 120% }\n	</style>\n</head>\n<body lang=\"fr-FR\" dir=\"ltr\" style=\"background: transparent\">\n<p align=\"center\" style=\"margin-bottom: 0cm; line-height: 100%\"><font size=\"4\" style=\"font-size: 16pt\"><b>Diversification\nHelp</b></font></p>\n<p align=\"justify\" style=\"margin-bottom: 0cm; line-height: 100%\"><br/>\n\n</p>\n<p align=\"justify\" style=\"margin-bottom: 0cm; line-height: 100%\">To\nuse <i>Diversification</i><span style=\"font-variant: normal\"><span style=\"font-style: normal\">,\nplease first open a tree file, next a file containing the\ncorresponding fossils table. </span></span><span style=\"font-variant: normal\"><span style=\"font-style: normal\">You\nmay need to</span></span><span style=\"font-variant: normal\"><span style=\"font-style: normal\">\nset</span></span><span style=\"font-variant: normal\"><span style=\"font-style: normal\">s</span></span><span style=\"font-variant: normal\"><span style=\"font-style: normal\">\n</span></span><span style=\"font-variant: normal\"><span style=\"font-style: normal\">the\nadditional</span></span><span style=\"font-variant: normal\"><span style=\"font-style: normal\">\nparameters, which are&nbsp;:</span></span></p>\n<ul>\n	<li/>\n<p align=\"justify\" style=\"margin-bottom: 0cm; font-variant: normal; font-style: normal; line-height: 100%\">\n	the <i>Origin</i>, i.e. the starting time of evolution,</p>\n	<li/>\n<p align=\"justify\" style=\"margin-bottom: 0cm; font-variant: normal; font-style: normal; line-height: 100%\">\n	the <i>End</i>, i.e. the contemporary time,</p>\n	<li/>\n<p align=\"justify\" style=\"margin-bottom: 0cm; font-variant: normal; font-style: normal; line-height: 100%\">\n	the number of samplings of fossil times, which are uniformly drawn\n	from the intervals provided in the fossils file.</p>\n</ul>\n<p align=\"justify\" style=\"margin-bottom: 0cm; font-variant: normal; font-style: normal; line-height: 100%\">\n<br/>\n\n</p>\n<p align=\"justify\" style=\"margin-bottom: 0cm; font-variant: normal; font-style: normal; line-height: 100%\">\nThe estimation of the diversification and fossil find rates is\nlaunched by the item <i>Compute</i> of the <i>Diversification</i>\nmenu or alternatively the corresponding item of the toolbar. The\ncomputation may be stopped by clicking on the <i>Stop</i> button\n(results displayed afterwards will be based only on the iterations\nperformed so far). The mean rates and their standard deviations among\nthe fossil times samplings are displayed as soon as the computation\nends.</p>\n</body>\n</html>"

MolTiWindow::MolTiWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MolTiWindow)
{
    ui->setupUi(this);
    graphModel = new QTableGraphModel(this);
    partModel = NULL;
    proxyModel = new QSortFilterProxyModel;
    ui->tableView_graphs->setModel(graphModel);
    ui->tableView_graphs->setItemDelegateForColumn(0, new QItemBigStringDelegate());
    ui->tableView_graphs->setItemDelegateForColumn(1, new QItemIntDelegate());
    ui->tableView_graphs->setItemDelegateForColumn(2, new QItemDoubleDelegate());
    ui->tableView_graphs->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    ui->treeView_communities->setModel(proxyModel);
    ui->treeView_communities->setItemDelegateForColumn(0, new QItemClassDelegate());
    ui->treeView_communities->setItemDelegateForColumn(1, new QItemIntDelegate());
    ui->treeView_communities->setItemDelegateForColumn(2, new QItemDoubleDelegate());
     helpDialog = NULL;
     multi = NULL;
    annot = NULL;
    sap = NULL;
    sig = NULL;
    info = NULL;
    dict = NULL;
    dictName = NULL;
    part = NULL;
}

MolTiWindow::~MolTiWindow()
{
    delete graphModel;
    if(multi != NULL)
        freeMultiGraph(multi);
    if(annot != NULL)
        freeAnnotation(annot);
    if(sap != NULL)
        freeStatAnnotationPartition(sap);
    if(sig != NULL) {
        free((void*)sig->start);
        free((void*)sig->table);
        free((void*)sig);
    }
    if(info != NULL)
        freeOntologyInfo(info);
    if(dict != NULL)
        freeLexiTree(dict);
    if(dictName != NULL)
        freeLexiTree(dictName);
    if(part != NULL)
        freePartition(part);
    delete ui;
}



void MolTiWindow::openAnnotation()
{
    QString fileName = QFileDialog::getOpenFileName(this);
    if (!fileName.isEmpty()) {
        FILE *fi;
        char* name = strdup(qPrintable(fileName));
        if((fi = fopen(name, "r"))) {
            if(annot != NULL)
                freeAnnotation(annot);
            annot = readAnnotation(fi);
            fclose(fi);
            ui->lineEdit_annotation->setText(fileName);
        } else {
            QMessageBox::warning(this, QString(tr("IO error")), (QString(tr("Unable to read\n"))).append(fileName) );
        }
    }
}


void MolTiWindow::openDescription()
{
    QString fileName = QFileDialog::getOpenFileName(this);
    if (!fileName.isEmpty()) {
        FILE *fi;
        char* name = strdup(qPrintable(fileName));
        if((fi = fopen(name, "r"))) {
            if(info != NULL) {
                freeOntologyInfo(info);
                freeLexiTree(dict);
            }
            info = readOntologyInfo(fi);
            fclose(fi);
            dict = getDictFromTable(info->id, info->size);
            ui->lineEdit_description->setText(fileName);
        } else {
           QMessageBox::warning(this, QString(tr("IO error")), (QString(tr("Unable to read\n"))).append(fileName) );
        }
    }
}




void MolTiWindow::addGraph()
{
    QFileDialog dialog(this);
    dialog.setWindowModality(Qt::WindowModal);
    dialog.setAcceptMode(QFileDialog::AcceptOpen);
    dialog.setFileMode(QFileDialog::ExistingFiles);
    QStringList files;
    if (dialog.exec())
        files = dialog.selectedFiles();
    else
        return;
    if (files.size()>0) {
        FILE *fi;
        for(int i=0; i<files.size(); i++)
        {
            char* name = strdup(qPrintable(files.at(i)));
            if((fi = fopen(name, "r"))) {
                TypeGraph *graph = readGraph(fi);
                fclose(fi);
                if(graph != NULL) {
                    graphModel->addGraph(graph, files.at(i));
                    ui->tableView_graphs->update();
                }
            } else {
                QMessageBox::warning(this, QString(tr("IO error")), (QString(tr("Unable to read\n"))).append(files.at(i)));
            }
        }
    }
}

void MolTiWindow::removeGraphs()
{
    QModelIndexList selected = ui->tableView_graphs->selectionModel()->selectedRows();
    std::sort(selected.begin(), selected.end());
    for(int i=selected.size()-1; i>=0; i--)
        graphModel->removeGraph(selected.at(i).row());
}

bool MolTiWindow::saveAnnotated()
{
    QFileDialog dialog(this);
    dialog.setWindowModality(Qt::WindowModal);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    QStringList files;
    if (dialog.exec())
        files = dialog.selectedFiles();
    else
        return false;
    char* name = strdup(qPrintable(files.at(0)));
    FILE *f;
    if((f = fopen(name, "w"))) {
        fprintSignificantTableEmpiricalAnais(f, *(sig), sap->sizeG, sap->nameA, info);
        fclose(f);
    } else {
        QMessageBox::warning(this, QString(tr("IO error")),  (QString(tr("Unable to write\n"))).append(files.at(0)));
        return false;
     }
    return true;
}

bool MolTiWindow::savePartition()
{
    QFileDialog dialog(this);
    dialog.setWindowModality(Qt::WindowModal);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    QStringList files;
    if (dialog.exec())
        files = dialog.selectedFiles();
    else
        return false;
    char* name = strdup(qPrintable(files.at(0)));
    FILE *f;
    if((f = fopen(name, "w"))) {
        fprintPartitionStandard(f, part, multi->name);
        fclose(f);
    } else {
        QMessageBox::warning(this, QString(tr("IO error")),  (QString(tr("Unable to write\n"))).append(files.at(0)));
        return false;
     }
    return true;
}

void MolTiWindow::help()
{
   if(helpDialog == NULL) {
        helpDialog = new QHelpDialog(this);
        QFile file(":/HelpMolTi.html");
        if(file.open(QIODevice::ReadOnly)) {
            helpDialog->setContent(QTextStream(&file).readAll());
            file.close();
        } else {
            helpDialog->setContent(QString(HELP_PAGE));
        }
    }
    helpDialog->show();
}

void MolTiWindow::about()
{
   QMessageBox::about(this, tr("About MolTi"),
            tr("<b>MolTi</b> detect communities from multiplex networks, i.e. sets of graphs sharing the same set of vertices but with different edges categories"));
}

void MolTiWindow::computePartition ()
{
    TypeGraph **tmp;
    int size = 0;
    QList<QGraphItem*> *list = graphModel->getList();
    tmp = (TypeGraph**) malloc((list->size())*sizeof(TypeGraph*));
    for(int i=0; i<list->size(); i++)
        if(list->at(i)->present)
            tmp[size++] = list->at(i)->graph;
    if(multi != NULL)
        freeMultiGraph(multi);
    multi = toMultiGraph(tmp, size);
    free((void*)tmp);
    if(size<1)
        return;
    double gamma = ui->doubleSpinBox_resolution->value();
    TypePartition *p = (TypePartition*) malloc(sizeof(TypePartition));
    *p = getPartition(multi, LouvainType, &(gamma));
    QPartitionModel *pm;
    if(annot != NULL) {
        if(sap != NULL) {
            freeStatAnnotationPartition(sap);
        }
        sap = newStatAnnotationPartition(annot, p, multi->name);
        double threshold = ui->doubleSpinBox_threshold->value();
        if(sig != NULL) {
            free((void*)sig->start);
            free((void*)sig->table);
            free((void*)sig);
        }
        sig = (TypeSignificantTable*) malloc(sizeof(TypeSignificantTable));
        *(sig) = getSignificantTableBonferroni(sap, threshold);
        double *val = (double*) malloc(p->sizeAtom*sizeof(double));
        for(int c=0; c<p->sizeAtom; c++)
            if(sig->start[c+1]>sig->start[c])
                val[c] = sig->table[sig->start[c]].corrected;
            else
                val[c] = -1.;
        pm = new QPartitionValueModel(p, multi->name, val, this);
    } else {
        pm = new QPartitionModel(p, multi->name, this);
    }
    proxyModel->setSourceModel(pm);
    ui->treeView_communities->sortByColumn(0, Qt::AscendingOrder);
    for (int column = 0; column < pm->columnCount(); ++column)
        ui->treeView_communities->resizeColumnToContents(column);
    //ui->treeView_communities->setModel(pm);
    connect(ui->treeView_communities->selectionModel(), SIGNAL(currentRowChanged(QModelIndex,QModelIndex)), this, SLOT(updateTextBrowser(QModelIndex,QModelIndex)));
    if(partModel != NULL)
        delete partModel;
    if(part != NULL)
        freePartition(part);
    partModel = pm;
    part = p;
    QStringList listName;
    for(int n=0; n<multi->sizeGraph; n++)
        listName.append(QString(multi->name[n]));
    QCompleter *completer = new QCompleter(listName, this);
    completer->setCaseSensitivity(Qt::CaseInsensitive);
    ui->lineEdit->setCompleter(completer);
    if(dictName != NULL)
        freeLexiTree(dictName);
    dictName = getDictFromTable(multi->name, multi->sizeGraph);

}

void MolTiWindow::goTo() {
    qDebug() << "go to" << ui->lineEdit->text();
    char* name = strdup(qPrintable(ui->lineEdit->text()));
    int n;
    if((n = findWordLexi(name, dictName)) >= 0) {
        QModelIndex index = proxyModel->mapFromSource(partModel->index(part->atom[n], 0, QModelIndex()));
        ui->treeView_communities->selectionModel()->select(index, QItemSelectionModel::ClearAndSelect| QItemSelectionModel::Rows);
        ui->treeView_communities->selectionModel()->currentRowChanged(index, QModelIndex());
        ui->treeView_communities->scrollTo( index, QAbstractItemView::EnsureVisible);
    }
}

#define SIZE_BUFFER_INFO 100000
void MolTiWindow::updateTextBrowser(const QModelIndex & current, const QModelIndex & previous)
{
    char buffer[SIZE_BUFFER_INFO], *tmp;
    QModelIndex index = proxyModel->mapToSource(current);
    if(!index.isValid() || !partModel->isValidIndex(index))
        return;
    tmp = buffer;
    if(partModel->isAtom(index))
    {
        int a = partModel->toAtom(index);
        tmp += sprintf (tmp, "Community %d\n", a+1);
        if(sig != NULL) {
            int i;
            for(i=sig->start[a]; i<sig->start[a+1]; i++) {
                int n;
                tmp += sprintf (tmp, "%s\t%.3lE\t%.3lE\t(%d %d %d %d)", sap->nameA[sig->table[i].ontology], sig->table[i].fisher, sig->table[i].corrected, sig->table[i].effCO, sig->table[i].effC-sig->table[i].effCO, sig->table[i].effO-sig->table[i].effCO, sap->sizeG-sig->table[i].effO-sig->table[i].effC+sig->table[i].effCO);
                if(info != NULL && (n = findWordLexi(sap->nameA[sig->table[i].ontology], dict)) >= 0)
                    tmp += sprintf (tmp, "\t%s", info->name[n]);
                tmp += sprintf (tmp, "\n");
            }
        }

    } else {
        int g = partModel->toItem(index);
        tmp += sprintf (tmp, "%s\n", multi->name[g]);
        if(sap != NULL) {
             int i;
            for(i=sap->startG[g]; i<sap->startG[g+1]; i++) {
                int n;
                tmp += sprintf (tmp, "\t%s", sap->nameA[sap->itemG[i]]);
                if(info != NULL && (n = findWordLexi(sap->nameA[sap->itemG[i]], dict)) >= 0)
                    tmp += sprintf (tmp, "\t%s", info->name[n]);
                tmp += sprintf (tmp, "\n");
            }
        }
    }
    ui->textBrowser->setText(buffer);
}

/*
#define SIZE_BUFFER_INFO 100000
void part_selected (GtkTreeView *tree_view, GtkTreePath *path, GtkTreeViewColumn *column, gpointer user_data) {
    GtkTreeModel *model;
    GtkTreeIter   iter;
    GtkTextBuffer *text_buffer;
    gchar buffer[SIZE_BUFFER_INFO];

    text_buffer = (GtkTextBuffer*) user_data;
    model = gtk_tree_view_get_model(tree_view);
    buffer[0] = '\0';
    switch(gtk_tree_path_get_depth (path)) {
        case 1:
            if (gtk_tree_model_get_iter(model, &iter, path)) {
                gchar *name, *tmp;
                int c;
                gtk_tree_model_get(model, &iter, COLUMN_PART_IDENT, &name, -1);
                tmp = name +strlen("atom ");
                c = atoi(tmp);
                tmp = buffer;
                tmp += sprintf (tmp, "Class n %d\n", c);
                if(sig != NULL) {
                    int i;
                    for(i=sig->start[c]; i<sig->start[c+1]; i++) {
                        int n;
                        tmp += sprintf (tmp, "%s\t%.3lE\t%.3lE\t(%d %d %d %d)", sap->nameA[sig->table[i].ontology], sig->table[i].fisher, sig->table[i].corrected, sig->table[i].effCO, sig->table[i].effC-sig->table[i].effCO, sig->table[i].effO-sig->table[i].effCO, sap->sizeG-sig->table[i].effO-sig->table[i].effC+sig->table[i].effCO);
                        if(info != NULL && (n = findWordLexi(sap->nameA[sig->table[i].ontology], dict)) >= 0)
                            tmp += sprintf (tmp, "\t%s", info->name[n]);
                        tmp += sprintf (tmp, "\n");
                    }
                }
                g_free(name);
            }
            break;
        case 2:
            if (gtk_tree_model_get_iter(model, &iter, path)) {
                gchar *name, *tmp;
                gtk_tree_model_get(model, &iter, COLUMN_PART_IDENT, &name, -1);
                tmp = buffer;
                tmp += sprintf (tmp, "%s\n", name);
                if(sap != NULL) {
                    int g;
                    TypeLexiTree *dict;
                    dict = getDictFromTable(sap->nameG, sap->sizeG);
                    if((g = findWordLexi(name, dict)) >= 0) {
                        int i;
                        for(i=sap->startG[g]; i<sap->startG[g+1]; i++) {
                            int n;
                            tmp += sprintf (tmp, "\t%s", sap->nameA[sap->itemG[i]]);
                            if(info != NULL && (n = findWordLexi(sap->nameA[sap->itemG[i]], dict)) >= 0)
                                tmp += sprintf (tmp, "\t%s", info->name[n]);
                            tmp += sprintf (tmp, "\n");
                        }
                    }
                }
                g_free(name);
            }
            break;
        default:
            break;
    }
    gtk_text_buffer_set_text (text_buffer, buffer, strlen(buffer));
}
*/
