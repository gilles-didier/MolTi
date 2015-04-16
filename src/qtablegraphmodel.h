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




#ifndef QTableGraphModel_H
#define QTableGraphModel_H
#include <QAbstractTableModel>
#include "Graph.h"


class QGraphItem
{
public:
    QGraphItem(TypeGraph *g, QString n, bool p, double d) {
        graph = g;
        name = n;
        present = p;
        density = d;
    };
    ~QGraphItem(){};
    void setPresent(bool b) {
        present = b;
    };
    TypeGraph *graph;
    QString name;
    bool present;
    double density;
};

class QTableGraphModel : public QAbstractTableModel
{
    Q_OBJECT
public:
    QTableGraphModel(QObject *parent=0);
    ~QTableGraphModel();
    int rowCount(const QModelIndex &parent = QModelIndex()) const ;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    QModelIndex index(int row, int column, const QModelIndex & parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
//    bool setData(const QModelIndex & index, const QVariant & value, int role = Qt::EditRole);
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;
    Qt::ItemFlags flags(const QModelIndex &index) const;
    void setHeaderRow(int section, QString s) const;
    void addGraph(TypeGraph* graph, QString name);
    void removeGraph(int r);
    QList<QGraphItem*> *getList();
 private:
    QList<QGraphItem*> *list;
};

#endif // QTableGraphModel_H
