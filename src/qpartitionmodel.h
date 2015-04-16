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




#ifndef QPartitionModel_H
#define QPartitionModel_H
#include <QAbstractTableModel>
#include "Partition.h"


class QPartitionModel : public QAbstractTableModel
{
    Q_OBJECT
public:
    QPartitionModel(TypePartition *partition, char **n = 0, QObject *parent=0);
    ~QPartitionModel();
    int rowCount(const QModelIndex &parent = QModelIndex()) const ;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    QModelIndex index(int row, int column, const QModelIndex & parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;
    bool hasChildren(const QModelIndex & parent = QModelIndex()) const;
    Qt::ItemFlags flags(const QModelIndex &index) const;
    QModelIndex parent(const QModelIndex &index) const;
    void setHeaderRow(int section, QString s) const;
    bool isValidIndex(const QModelIndex &index) const;
    bool isAtom(const QModelIndex &index) const;
    bool isItem(const QModelIndex &index) const;
    int toAtom(const QModelIndex &index) const;
    int toItem(const QModelIndex &index) const;
private:
    TypePartition *part;
    TypePartitionCompact *pc;
    char **name;
};

#endif // QPartitionModel_H
