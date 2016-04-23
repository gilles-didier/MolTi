#include <QtWidgets>
#include "qhelpdialog.h"

QHelpDialog::QHelpDialog(QWidget * parent) :
    QDialog(parent, Qt::Dialog)
{

    setupUi(this);
}

QHelpDialog::~QHelpDialog()
{

}

void QHelpDialog::setContent(const QUrl &url) {
    textBrowser->setSource(url);
    textBrowser->show();
}

