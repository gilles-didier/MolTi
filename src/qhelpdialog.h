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




#ifndef QHELPDIALOG_H
#define QHELPDIALOG_H
#include <QDialog>
#include "ui_qhelpdialog.h"

class QHelpDialog : public QDialog, private Ui::MyHelpDialog
{
public:
    QHelpDialog(QWidget * parent = 0);
    ~QHelpDialog();
    void setContent(const QString &text);
    void setBaseUrl(const QUrl &url);
private:
    QUrl baseUrl;
};

#endif // QHELPDIALOG_H
