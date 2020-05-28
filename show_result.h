#ifndef SHOW_RESULT_H
#define SHOW_RESULT_H

#include <QDialog>
#include "data.h"

namespace Ui {
class Show_result;
}

class Show_result : public QDialog
{
    Q_OBJECT

public:
    explicit Show_result(QWidget *parent, Data* &d);
    ~Show_result();

private:
    Ui::Show_result *ui;
    Data* data{};
};

#endif // SHOW_RESULT_H
