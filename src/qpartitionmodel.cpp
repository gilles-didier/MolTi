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
#include "qpartitionmodel.h"
#include <QDebug>

#define NOTSET -1.
#define EPS 0.001

QPartitionModel::QPartitionModel(TypePartition *p, char **n, QObject *parent):QAbstractTableModel(parent)
{
    part = p;
    pc = getPartitionCompact(part);
    name = n;
}

QPartitionModel::~QPartitionModel()
{
    freePartitionCompact(pc);
}

int QPartitionModel::rowCount(const QModelIndex &parent) const
{
    if(parent.isValid()) {
        if(isAtom(parent))
            return pc->start[parent.internalId()+1]-pc->start[parent.internalId()];
        else
            return 0;
    } else
        return part->sizeAtom;
}

bool QPartitionModel::isAtom(const QModelIndex &index) const
{
    return index.isValid() && (((int)index.internalId())<part->sizeAtom);
}

bool QPartitionModel::isItem(const QModelIndex &index) const
{
    return index.isValid() && (((int)index.internalId())>=part->sizeAtom);
}

int QPartitionModel::toAtom(const QModelIndex &index) const
{
    if(index.isValid())
        return index.internalId();
    else return -1;
}

int QPartitionModel::toItem(const QModelIndex &index) const
{
    if(index.isValid() && isItem(index))
        return index.internalId()-part->sizeAtom;
    else return -1;
}
int QPartitionModel::columnCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return 1;
    else
        return 2;
}
bool QPartitionModel::isValidIndex(const QModelIndex &index) const {
    return index.internalId()>=0 && index.internalId()<(part->sizeAtom+part->sizeItem);
}

Qt::ItemFlags QPartitionModel::flags(const QModelIndex &index) const
{
    if (!index.isValid())
        return 0;
    if(isAtom(index))
        return Qt::ItemIsSelectable | Qt::ItemIsEnabled;
    else
        return Qt::ItemIsSelectable | Qt::ItemIsEnabled;
}

QModelIndex QPartitionModel::index(int row, int column, const QModelIndex & parent) const
{
    if (parent.isValid()) {
        return createIndex(row, column, part->sizeAtom+pc->item[pc->start[parent.internalId()]+row]);
    } else {
        return createIndex(row, column, row);
    }

}

QModelIndex QPartitionModel::parent(const QModelIndex &index) const
{
    if (!index.isValid())
        return QModelIndex();
    if(isAtom(index))
        return QModelIndex();
    int a = part->atom[index.internalId()-part->sizeAtom];
    return createIndex(a, 0, a);
}

bool QPartitionModel::hasChildren(const QModelIndex & index) const
{
    if (!index.isValid())
        return true;
    return isAtom(index);
}


QVariant QPartitionModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();
    if(((int)index.internalId())<part->sizeAtom) {
        int a = index.internalId();
        switch(role){
            case Qt::DisplayRole:
                switch(index.column()) {
                    case 0: return QVariant((int)a+1);
                    case 1: return QVariant(pc->start[a+1]-pc->start[a]);
                    default: return QVariant();
                }
            case Qt::TextAlignmentRole:
                switch(index.column()) {
                    case 0: return Qt::AlignLeft + Qt::AlignVCenter;
                    case 1: return Qt::AlignRight + Qt::AlignVCenter;
                    default: return QVariant();
                }
        }
    } else {
        int g = index.internalId()-part->sizeAtom;
        switch(role){
            case Qt::DisplayRole:
                switch(index.column()) {
                    case 0:
                        if(name != NULL && name[g] != NULL)
                            return QString(name[g]);
                        else
                            return QString("gene ").append(QString::number(g));
                     default: return QVariant();
                }
            case Qt::TextAlignmentRole:
                switch(index.column()) {
                    case 0: return Qt::AlignLeft + Qt::AlignVCenter;
                    default: return QVariant();
                }
        }
    }
    return QVariant();
}

QVariant QPartitionModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    switch(orientation) {
        case Qt::Horizontal:
            switch(role){
                case Qt::DisplayRole:
                    switch(section)
                    {
                        case 0: return QString("Class");
                        case 1: return QString("Size");
                        default: return QVariant();
                    }
                    break;
                default: return QVariant();
            }
        default: return QVariant();
    }
}
