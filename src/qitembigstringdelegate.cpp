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




#include "qitembigstringdelegate.h"
#include <QPainter>
#include <QFontMetrics>

QItemBigStringDelegate::QItemBigStringDelegate() : QStyledItemDelegate()
{
    font= QFont("Times", -1, QFont::Normal);
    alpha = 0.4;
}

QItemBigStringDelegate::~QItemBigStringDelegate()
{

}

void QItemBigStringDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option,
                         const QModelIndex &index) const
{
    if (index.data().canConvert<QString>()) {
        QString st = index.data().toString();
        QString news = QString(st);
        int w = painter->fontMetrics().width(news);
        if(w>option.rect.width()) {
             int p;
             bool bs = true;
             p = alpha*news.length();
             news.insert(p,QString("..."));
             while(painter->fontMetrics().width(news)>option.rect.width() && (p>0 || news.length()>3)) {
                 if(bs && p>0) {
                    news.remove(--p,1);
                    bs = false;
                } else {
                     news.remove(p+3,1);
                     bs = true;
                }
            }
        }
        painter->drawText(option.rect, qvariant_cast<int>(index.data(Qt::TextAlignmentRole)), news);
    } else {
        QStyledItemDelegate::paint(painter, option, index);
    }
}
