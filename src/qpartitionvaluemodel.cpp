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




#include "qpartitionvaluemodel.h"
#include "float.h"
#include <QLocale>

QPartitionValueModel::QPartitionValueModel(TypePartition *p, char **n, double *v, QObject *parent)
    :QPartitionModel(p, n, parent)
{
    val = v;
}

QPartitionValueModel::~QPartitionValueModel()
{

}

int QPartitionValueModel::columnCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return 1;
    else
        return (val != NULL)?3:2;
}

QVariant QPartitionValueModel::data(const QModelIndex &index, int role) const
{
    QLocale locale; // Constructs a default QLocale
    if(index.column()<2)
        return QPartitionModel::data(index, role);
    if (!index.isValid() || isItem(index) || val == NULL)
            return QVariant();
    switch(role) {
        case Qt::DisplayRole:
            if(val[index.internalId()]>=0)
//                return locale.toString(val[index.internalId()], 'e', 2);
                return val[index.internalId()];
            else
                return DBL_MAX;
        case Qt::TextAlignmentRole:
            return Qt::AlignRight + Qt::AlignVCenter;
    }
    return QVariant();
}

QVariant QPartitionValueModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    switch(orientation) {
        case Qt::Horizontal:
            switch(role){
                case Qt::DisplayRole:
                    switch(section)
                    {
                        case 0: return QString("Class");
                        case 1: return QString("Size");
                        case 2: return QString("q-value");
                        default: return QVariant();
                    }
                    break;
                default: return QVariant();
            }
        default: return QVariant();
    }
}
