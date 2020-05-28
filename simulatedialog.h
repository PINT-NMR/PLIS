#ifndef SIMULATEDIALOG_H
#define SIMULATEDIALOG_H

#include "data.h"
#include <QVector>
#include <QDialog>
#include "calculations.h"


namespace Ui {
class SimulateDialog;
}

class SimulateDialog : public QDialog
{
    Q_OBJECT

public:
    explicit SimulateDialog(QWidget *parent = nullptr);
    explicit SimulateDialog(QWidget *parent, QString mode);
    ~SimulateDialog();

    QVector<Data*> dataSets_local;
    QString mode;
    int current_dataSet{};
    void setup(QVector<Data *> &datasets, QString mode);

private slots:

    void on_simulateButton_clicked();

    void on_responce_s_valueChanged(double arg1);

    void on_responce_I_valueChanged(double arg1);

    void on_responce_L_valueChanged(double arg1);

    void on_responce_start_valueChanged(double arg1);

    void on_responce_end_valueChanged(double arg1);

    void on_cancelButton_clicked();

private:
    Ui::SimulateDialog *ui;
    QDoubleValidator *v=new QDoubleValidator(0.001,99999,3,this);

    void simulateGenerateData(QVector<double> &x, double noise);

    void one_site_sim();
    void two_site_sim();
    void comp_sim();
};

#endif // SIMULATEDIALOG_H
