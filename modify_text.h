#ifndef MODIFY_TEXT_H
#define MODIFY_TEXT_H

#include <QDialog>
#include "qcustomplot.h"

namespace Ui {
class Modify_text;
}

class Modify_text : public QDialog
{
    Q_OBJECT

public:
    Modify_text(QWidget *parent, QCPAxis* axis, QFont &x, QFont &y, QString &xlabel, QString &ylabel);
    Modify_text(QWidget *parent, QCPTextElement* &title);
    ~Modify_text();

private slots:
    void on_font_clicked();

    void on_label_text_clicked();

    void on_both_clicked();

private:
    Ui::Modify_text *ui;
    QCPAxis *axis_local{};
    QFont *xlabel_font{},*ylabel_font{};
    QString *xlabel{}, *ylabel{};
    QCPTextElement* title_local{};
    int mode{};
};

#endif // MODIFY_TEXT_H
