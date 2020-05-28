#ifndef GUESS_PARAMETERS_H
#define GUESS_PARAMETERS_H

#include <QDialog>
#include "data.h"
#include "calculations.h"

namespace Ui {
class Guess_parameters;
}

class Guess_parameters : public QDialog
{
    Q_OBJECT

public:
    explicit Guess_parameters(QWidget *parent, Data *&d);
    ~Guess_parameters();

    bool a0{},a1{},a2{},a3{},a4{},a5{};

private slots:
    void on_pushButton_clicked();

    void on_kd1_check_one_stateChanged(int arg1);

    void on_start_check_one_stateChanged(int arg1);

    void on_end_check_one_stateChanged(int arg1);

    void on_kd1_check_two_stateChanged(int arg1);

    void on_kd2_check_two_stateChanged(int arg1);

    void on_start_check_two_stateChanged(int arg1);

    void on_intermediate_check_two_stateChanged(int arg1);

    void on_end_check_two_stateChanged(int arg1);

    void on_kd_min_valueChanged(double arg1);

    void on_kd_max_valueChanged(double arg1);

    void on_start_min_valueChanged(double arg1);

    void on_start_max_valueChanged(double arg1);

    void on_end_min_valueChanged(double arg1);

    void on_end_max_valueChanged(double arg1);

    void on_kd_check_stateChanged(int arg1);

    void on_start_check_range_stateChanged(int arg1);

    void on_end_check_range_stateChanged(int arg1);

    void on_kd_min_2_valueChanged(double arg1);

    void on_kd_max_2_valueChanged(double arg1);

    void on_kd2_min_valueChanged(double arg1);

    void on_kd2_max_valueChanged(double arg1);

    void on_start_min_2_valueChanged(double arg1);

    void on_start_max_2_valueChanged(double arg1);

    void on_intermediate_min_valueChanged(double arg1);

    void on_intermediate_max_valueChanged(double arg1);

    void on_end_min_2_valueChanged(double arg1);

    void on_end_max_2_valueChanged(double arg1);

    void on_kd_check_range_2_stateChanged(int arg1);

    void on_kd2_check_range_stateChanged(int arg1);

    void on_start_check_range_2_stateChanged(int arg1);

    void on_intermediate_check_range_stateChanged(int arg1);

    void on_end_check_range_2_stateChanged(int arg1);

    void on_cancel_clicked();

    void on_kd1_check_two_2_stateChanged(int arg1);

    void on_kd2_check_two_2_stateChanged(int arg1);

    void on_kdc_check_stateChanged(int arg1);

    void on_start_check_two_2_stateChanged(int arg1);

    void on_end_check_two_2_stateChanged(int arg1);

    void on_kd_check_range_comp_stateChanged(int arg1);

    void on_kd2_check_range_comp_stateChanged(int arg1);

    void on_kdc_check_range_stateChanged(int arg1);

    void on_start_check_range_comp_stateChanged(int arg1);

    void on_end_check_range_comp_stateChanged(int arg1);

    void on_kd_min_comp_valueChanged(double arg1);

    void on_kd_max_comp_valueChanged(double arg1);

    void on_kd2_min_comp_valueChanged(double arg1);

    void on_kd2_max_comp_valueChanged(double arg1);

    void on_kdc_min_valueChanged(double arg1);

    void on_kdc_max_valueChanged(double arg1);

    void on_start_min_comp_valueChanged(double arg1);

    void on_start_max_comp_valueChanged(double arg1);

    void on_end_min_comp_valueChanged(double arg1);

    void on_end_max_comp_valueChanged(double arg1);

    void on_r20_check_stateChanged(int arg1);

    void on_kex_check_stateChanged(int arg1);

    void on_pb_check_stateChanged(int arg1);

    void on_dw_check_stateChanged(int arg1);

    void on_r20_min_valueChanged(double arg1);

    void on_r20_max_valueChanged(double arg1);

    void on_kex_min_valueChanged(double arg1);

    void on_kex_max_valueChanged(double arg1);

    void on_pb_min_valueChanged(double arg1);

    void on_pb_max_valueChanged(double arg1);

    void on_dw_min_valueChanged(double arg1);

    void on_dw_max_valueChanged(double arg1);

    void on_dw_check_r_stateChanged(int arg1);

    void on_kex_check_r_stateChanged(int arg1);

    void on_pb_check_r_stateChanged(int arg1);

    void on_r20_check_r_stateChanged(int arg1);

    void on_tabWidget_one_currentChanged(int index);

    void on_tabWidget_two_currentChanged(int index);

    void on_tabWidget_comp_currentChanged(int index);

    void on_tabWidget_cpmg_currentChanged(int index);

private:
    Ui::Guess_parameters *ui;
    Data *data{};
};

#endif // GUESS_PARAMETERS_H
