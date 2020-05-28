#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <QDialog>
#include "data.h"
#include "fitmrq.h"
#include "jackknife.h"

namespace Ui {
class Optimize;
}

class Optimize : public QDialog
{
    Q_OBJECT

public:
    explicit Optimize(QWidget *parent, Data* &d);
    ~Optimize();

private slots:
    void on_pushButton_3_clicked();

    void on_apply_clicked();

    void on_pushButton_2_clicked();

private:
    Ui::Optimize *ui;
    Data * data{};
    Data * last=nullptr;
};

#endif // OPTIMIZE_H
