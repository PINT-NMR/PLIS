#ifndef MODIFY_DATA_H
#define MODIFY_DATA_H

#include <QDialog>
#include "qcustomplot.h"
#include "data.h"

namespace Ui {
class Modify_data;
}

class Modify_data : public QDialog
{
    Q_OBJECT

public:
    Modify_data(QWidget *parent, QVector<QCPScatterStyle::ScatterShape> s, Data *&d);
    ~Modify_data();

    int selected_shape{};
    QPen pen{};
    QPen pen_curve{};

private slots:

    void on_okButton_clicked();

    void on_cancelButton_clicked();

private:
    Ui::Modify_data *ui;
    QVector<QCPScatterStyle::ScatterShape> shapes;
    Data* data{};

};

#endif // MODIFY_DATA_H
