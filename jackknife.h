#ifndef JACKKNIFE_H
#define JACKKNIFE_H

#include "data.h"
#include "fitmrq.h"
#include "calculations.h"


class Jackknife
{
public:
    Jackknife(Data* &data, const std::vector<bool> &_ia);

    Jackknife(Data* &data, const std::vector<bool> &_ia, QString _mode);

    void compute();

    void compute_cpmg();

    std::vector<std::vector<double>> a{};
    double deviation2{},deviation_kd1{},deviation_kd2{};

private:
    void calculate();
    void clearVectors();
    void writeJackDataSet(int i);
    std::vector<double> calcAverageA(double ndata);
    std::vector<double> calcError();
    void assignErrorsToData(const std::vector<double> &jackError);

    Data* dataset;
    std::vector<bool> ia{};
    QString mode{"standard"};
    QVector<double> responce{};
    QVector<double> conc{}, protein{}, error{}, volume{}, comp{},r2eff{},n_cpmg{};

};

#endif // JACKKNIFE_H
