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




#ifndef QPARTITIONITEMDELEGATE_H
#define QPARTITIONITEMDELEGATE_H
#include <QStyledItemDelegate>
#include <QString>


class QPartitionItemDelegate : public QStyledItemDelegate
{
public:
    QPartitionItemDelegate();
    ~QPartitionItemDelegate();
    QString displayText(const QVariant & value, const QLocale & locale) const;
};

#endif // QPARTITIONITEMDELEGATE_H
