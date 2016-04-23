#ifndef QHELPDIALOG_H
#define QHELPDIALOG_H
#include <QDialog>
#include "ui_qhelpdialog.h"

class QHelpDialog : public QDialog, private Ui::MyHelpDialog
{
public:
    QHelpDialog(QWidget * parent = 0);
    ~QHelpDialog();
    void setContent(const QUrl &url);
    void setBaseUrl(const QUrl &url);
};

#endif // QHELPDIALOG_H
