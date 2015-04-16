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




#ifndef MOLTIWINDOW_H
#define MOLTIWINDOW_H

#include <QString>
#include <QMainWindow>
#include <QSortFilterProxyModel>
#include "Utils.h"
#include "Graph.h"
#include "Partition.h"
#include "Louvain.h"
#include "EdgesComposition.h"
#include "Annotation.h"
#include "StatAnnotationPartition.h"
#include "qtablegraphmodel.h"
#include "qpartitionmodel.h"
#include "qpartitionvaluemodel.h"
#include "qhelpdialog.h"

namespace Ui {
class MolTiWindow;
}

class MolTiWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MolTiWindow(QWidget *parent = 0);
    ~MolTiWindow();
private slots:
    void openAnnotation();
    void openDescription();
    void addGraph();
    void removeGraphs();
    void computePartition();
    bool saveAnnotated();
    bool savePartition();
    void help();
    void about();
    void updateTextBrowser(const QModelIndex & current, const QModelIndex & previous);
    void goTo();
private:
    Ui::MolTiWindow *ui;
    TypeMultiGraph *multi;
    TypeAnnotation *annot;
    TypeStatAnnotationPartition *sap;
    TypeSignificantTable *sig;
    TypeOntologyInfo *info;
    TypeLexiTree *dict;
    TypePartition *part;
    QTableGraphModel *graphModel;
    QPartitionModel *partModel;
    QSortFilterProxyModel *proxyModel;
    QHelpDialog *helpDialog;
    TypeLexiTree *dictName;
};

#endif // MOLTIWINDOW_H
