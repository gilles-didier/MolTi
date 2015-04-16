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




#include "qpartitionitemdelegate.h"

QPartitionItemDelegate::QPartitionItemDelegate()
{

}

QPartitionItemDelegate::~QPartitionItemDelegate()
{

}


QString QPartitionItemDelegate::displayText(const QVariant & value, const QLocale & locale) const
{
    bool ok;
    double v = value.toDouble(&ok);
//    QLocale locale; // Constructs a default QLocale
    if(ok)
        return locale.toString(v, 'e', 2);
    else
        return value.toString();
}
