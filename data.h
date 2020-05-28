#ifndef DATA_H
#define DATA_H

#include <string>
#include <QVector>
#include "qcustomplot.h"

//namespace std {
class Data;
//}
class Data
{
public:
    Data();
    ~Data()=default;
    void read_Data(std::string &);
    void change_unit_of_data(QString &arg1, QString data);
    void update_protein_conc(int row=-1);
    void update_comp_vector(double start);
    void change_unit_of_responce(QString);

    QVector<double> responceVector{};
    QVector<double> concVector{};
    QVector<double> volumeVector{};
    QVector<double> n_cpmgVector{}, R2effVector{},dyVector{};

    std::vector<double> a{0,0,0}, dyda{0,0,0};
    QVector<QVector<double>> calc_data{QVector<double>{},QVector<double>{}};
    QVector<double> protein_conc_vector{}, comp_vector{}, error_vector{};
    double protein_conc{};
    double ligand_cpmg{}, kd_cpmg{};
    double kd_error{}, kd1_error{}, kd2_error{}, kdc_error{};
    double chi2{};
    QString ligand_unit{};
    QString unit{};
    QString name{};
    QString model{};

    bool data_visible;
    bool curve_visible;
    bool has_curve;
    int num_bind_site{};
    QCPScatterStyle style{};
    QPen pen{};
    QPen pen_curve{};

private:
    void change_data(double, QString);

};

#endif // DATA_H
