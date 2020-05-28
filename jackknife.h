#ifndef JACKKNIFE_H
#define JACKKNIFE_H

#include "data.h"
#include "fitmrq.h"
#include "calculations.h"


class Jackknife
{
public:
    Jackknife(Data* &data);
    void compute();

    void compute_cpmg();

    std::vector<std::vector<double>> a{std::vector<double>{}, std::vector<double>{}, std::vector<double>{}, std::vector<double>{}, std::vector<double>{}};
    double deviation2{},deviation_kd1{},deviation_kd2{};


private:
    void calculate();

    Data* dataset;
    QVector<double> responce{};
    QVector<double> conc{}, protein{}, error{}, volume{}, comp{},r2eff{},n_cpmg{};

};

#endif // JACKKNIFE_H
