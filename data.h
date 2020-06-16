#ifndef DATA_H
#define DATA_H

// 200614
// - changed meaning of dyda

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
    void clearVectors();
    void pushBackData(double response, double conc, double volume, double protConc, double err);
    void insertData(int i, double response, double conc, double volume, double protConc, double err);
    void removeData(int i);
    void pushBackCPMGdata(double n_cpmg, double r2eff, double dr2eff);
    void insertCPMGdata(int i, double n_cpmg, double r2eff, double dr2eff);
    void removeCPMGdata(int i);

    void addCPMGresults(const std::vector<double> &fitted_a, double fitted_chi2, const std::vector<std::vector<double>> &temp_vector);
    void removeFittedCurve();

    void setModelIfNoFit();
    void setModel_a_dydaOneSite();
    void setModel_a_dydaTwoSite();
    void setModel_a_dydaFourSite();
    void setModel_a_dydaCompTwoSite(double comp_conc);
    void setModel_a_dydaCPMG();

    void setResult(const std::vector<double> &_a, double _chi2, std::vector<double> &_calc_data0, std::vector<double> &_calc_data1);

    void dataAndCurveVisible(bool dataVis, bool curveVis);

    void setPenColor(int i);

    QVector<double> responceVector{};
    QVector<double> concVector{};
    QVector<double> volumeVector{};
    QVector<double> n_cpmgVector{}, R2effVector{},dyVector{};

    // N.B. dyda is not necessary so I will fill it with da (errors in a) instead
    // should rename later
    std::vector<double> a{0,0,0}, dyda{0,0,0};
    QVector<QVector<double>> calc_data{QVector<double>{},QVector<double>{}};
    QVector<double> protein_conc_vector{}, comp_vector{}, error_vector{};
    double protein_conc{};
    double ligand_cpmg{}, kd_cpmg{};
    double kd_error{}, kd1_error{}, kd2_error{}, kd3_error{}, kd4_error{}, kdc_error{};
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
