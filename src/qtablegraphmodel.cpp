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




#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cfloat>
#include "qtablegraphmodel.h"
#include <QDebug>

#define NOTSET -1.
#define EPS 0.001

QTableGraphModel::QTableGraphModel(QObject *parent):QAbstractTableModel(parent)
{
    list = new QList<QGraphItem*>();
}

QTableGraphModel::~QTableGraphModel()
{
    delete list;
}


int QTableGraphModel::rowCount(const QModelIndex & /*parent*/) const
{
    if(list != 0)
        return list->size();
    else
        return 0;
}

int QTableGraphModel::columnCount(const QModelIndex & /*parent*/) const
{
    return 3;
}

QModelIndex QTableGraphModel::index(int row, int column, const QModelIndex & parent) const
{
   return createIndex(row, column, (quintptr) row);
}

Qt::ItemFlags QTableGraphModel::flags(const QModelIndex &index) const
{
    if (!index.isValid())
        return 0;
    if(index.column() == 0)
        return Qt::ItemIsSelectable | Qt::ItemIsEnabled;// | Qt::ItemIsUserCheckable | Qt::ItemIsTristate | Qt::ItemIsEditable;// | Qt::ItemIsEditable;
    else
        return Qt::ItemIsSelectable | Qt::ItemIsEnabled;
}

void QTableGraphModel::addGraph(TypeGraph* graph, QString name)
{
    beginInsertRows(QModelIndex(), rowCount(), rowCount());
    list->append(new QGraphItem(graph, name, true, getGraphDensity(graph)));
    endInsertRows();
}

void QTableGraphModel::removeGraph(int r)
{
    beginRemoveRows(QModelIndex(), r, r);
    freeGraph(list->at(r)->graph);
    list->removeAt(r);
    endRemoveRows();
}

QList<QGraphItem*> *QTableGraphModel::getList() {
    return list;
}

QVariant QTableGraphModel::data(const QModelIndex &index, int role) const
{
    switch(role){
        case Qt::DisplayRole:
            switch(index.column())
            {
                case 0: return list->at(index.row())->name;
                case 1: return list->at(index.row())->graph->sizeGraph;
                case 2: return list->at(index.row())->density;
                default: return QVariant();
            }
            break;
  /*       case Qt::CheckStateRole:
            switch(index.column())
            {
                case 0: return QVariant(list->at(index.internalId())->present);//Qt::Checked;
                default: return QVariant();
            }
            break;
       case Qt::EditRole:
            switch(index.column())
            {
                case 0: return list->at(index.internalId()).present;
                default: return QVariant();
            }
            break;
 */       case Qt::TextAlignmentRole:
            switch(index.column())
            {
                case 0: return Qt::AlignLeft + Qt::AlignVCenter;
                case 1: return Qt::AlignRight + Qt::AlignVCenter;
                case 2: return Qt::AlignRight + Qt::AlignVCenter;
                default: return Qt::AlignRight + Qt::AlignVCenter;
            }
    }
    return QVariant();
}

/*
bool QTableGraphModel::setData(const QModelIndex & index, const QVariant & value, int role) {
    if (!index.isValid())
        return true;
    int g = (int) index.internalId();
    if (index.column()==0 && role == Qt::CheckStateRole)
    {
//        list->at(g).present = (Qt::CheckState)value.toBool();
       list->operator[](g)->present = (Qt::CheckState)value.toBool();
    emit dataChanged(index,index);
       return true;
    }
    return QAbstractItemModel::setData ( index,value,role );
}
*/

QVariant QTableGraphModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    switch(orientation) {
        case Qt::Horizontal:
            switch(role){
                case Qt::DisplayRole:
                    switch(section)
                    {
                        case 0: return QString("File name");
                        case 1: return QString("# of vertices");
                        case 2: return QString("Density");
                        default: return QVariant();
                    }
                    break;
                default: return QVariant();
            }
        default: return QVariant();
    }
}
