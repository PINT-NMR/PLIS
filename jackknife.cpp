#include "jackknife.h"

// 200530
// - the ia vector from fitmrq is now in the constructor
//   and correct values of it are used to calculate the errors
// - some cosmetic changes

Jackknife::Jackknife(Data* &data, const std::vector<bool> &_ia) : dataset(data), ia(_ia)
{
    for (int i=0; i<signed(dataset->a.size()); ++i)
        a.push_back({});
}

Jackknife::Jackknife(Data* &data, const std::vector<bool> &_ia, QString _mode) :
    dataset(data), ia(_ia), mode(_mode)
{
    for (int i=0; i<signed(dataset->a.size()); ++i)
        a.push_back({});
    std::cerr << "size " << a.size() << "\n";
}


void Jackknife::compute() {
    for (int i=0;i<dataset->responceVector.size();i++) {
        writeJackDataSet(i);
        calculate();
    }
    std::vector<double> jackError = calcError();
    assignErrorsToData(jackError);
}

void Jackknife::clearVectors() {
    responce.clear();
    conc.clear();
    protein.clear();
    error.clear();
    volume.clear();
    comp.clear();
}

void Jackknife::writeJackDataSet(int i) {
    clearVectors();
    for (int j=0; j<dataset->responceVector.size(); j++)
        if (i!=j) {
            responce.push_back(dataset->responceVector[j]);
            conc.push_back(dataset->concVector[j]);
            protein.push_back(dataset->protein_conc_vector[j]);
            error.push_back(dataset->error_vector[j]);
            volume.push_back(dataset->volumeVector[j]);
            if(dataset->num_bind_site==3)
                comp.push_back(dataset->comp_vector[j]);
        }
}

std::vector<double> Jackknife::calcAverageA(double ndata) {
    std::vector<double> average(signed (a.size()),0.);
    for (unsigned int j=0;j<unsigned (dataset->a.size());j++){
        for (unsigned int k=0;k<ndata;k++)
            average[signed (j)]+=a[j][k];
        average[signed (j)] /= ndata;
    }
    return average;
}

std::vector<double> Jackknife::calcError() {
    double njack=a[0].size()-1;
    double ndata=a[0].size();
    std::vector<double> average = calcAverageA(ndata);
    std::vector<double> sumvec(a.size(), 0.);
    for (unsigned int j=0; j<unsigned (dataset->a.size()); ++j) {
        for (unsigned int i=0;i<ndata; i++)
            sumvec[j]+=(a[j][i]-average[j])*(a[j][i]-average[j]);
        sumvec[j] *= (njack/ndata);
    }
    return sumvec;
}

void Jackknife::assignErrorsToData(const std::vector<double> &jackError) {
    dataset->dyda = jackError;
    if (dataset->num_bind_site==1)
        dataset->kd_error=sqrt(jackError[2]);
    else if (dataset->num_bind_site==2) {
        dataset->kd1_error=sqrt(jackError[3]);
        dataset->kd2_error=sqrt(jackError[4]);
    }
    else if (dataset->num_bind_site==3) {
        dataset->kd1_error=sqrt(jackError[2]);
        dataset->kd2_error=sqrt(jackError[3]);
        dataset->kdc_error=sqrt(jackError[4]);
    }
    else if (dataset->num_bind_site==5) {
        dataset->kd1_error=sqrt(jackError[5]);
        dataset->kd2_error=sqrt(jackError[6]);
        dataset->kd3_error=sqrt(jackError[7]);
        dataset->kd4_error=sqrt(jackError[8]);
    }
}

void Jackknife::calculate() {
    std::vector<double> temp_a{dataset->a};
    if (dataset->num_bind_site==1) {
        void (*voidptr)(const double, const double, const std::vector<double> &, double &,
                        std::vector<double> &){fit_one_site_dilution};
        double (*dblptr)(const double, const double, const std::vector<double> &){fit_one_site_dilution};
        if (mode=="chemshift") {
            voidptr = fit_one_site;
            dblptr = fit_one_site;
        }
        Fitmrq calc(conc.toStdVector(),protein.toStdVector(),responce.toStdVector(),
                    volume.toStdVector(),error.toStdVector(),temp_a, ia, voidptr,dblptr);
        calc.fit();
        for (unsigned int i=0; i<a.size(); ++i)
            a[i].push_back(calc.a[i]);
    }
    else if(dataset->num_bind_site==2) {
        void (*voidptr)(const double, const double, const std::vector<double> &, double &,
                        std::vector<double> &){fit_two_sites_dilution};
        double (*dblptr)(const double, const double, const std::vector<double> &){fit_two_sites_dilution};
        if (mode=="chemshift") {
            voidptr = fit_two_sites;
            dblptr = fit_two_sites;
        }
        Fitmrq calc(conc.toStdVector(),protein.toStdVector(),responce.toStdVector(),
                    volume.toStdVector(),error.toStdVector(),temp_a,ia, voidptr, dblptr);;
        calc.fit();
        for (unsigned int i=0; i<a.size(); ++i) a[i].push_back(calc.a[i]);
    }
    else if (dataset->num_bind_site==3) {
        void (*voidptr)(const double, const double, const double, const std::vector<double> &, double &,
                        std::vector<double> &){fit_comp_dilution};
        double (*dblptr)(const double, const double, const double, const std::vector<double> &){fit_comp_dilution};
        if (mode=="chemshift") {
            voidptr = fit_comp;
            dblptr = fit_comp;
        }
        Fitmrq2 calc(conc.toStdVector(),protein.toStdVector(),comp.toStdVector(),responce.toStdVector(),
                    volume.toStdVector(),error.toStdVector(), temp_a,ia, voidptr, dblptr);
        calc.fit();
        for (unsigned int i=0; i<a.size(); ++i) a[i].push_back(calc.a[i]);
    }
    else if (dataset->num_bind_site==5) {
        void (*voidptr)(const double, const double, const std::vector<double> &, double &,
                        std::vector<double> &){fit_four_sites_dilution};
        double (*dblptr)(const double, const double, const std::vector<double> &){fit_four_sites_dilution};
        if (mode=="chemshift") {
            voidptr = fit_four_sites;
            dblptr = fit_four_sites;
        }
        Fitmrq calc(conc.toStdVector(), protein.toStdVector(), responce.toStdVector(),
                    volume.toStdVector(), error.toStdVector(), temp_a, ia, voidptr, dblptr);
        calc.fit();
        for (unsigned int i=0; i<a.size(); ++i) a[i].push_back(calc.a[i]);
    }
    else if (dataset->num_bind_site==4) { //CPMG
        Fitmrq3 calc(n_cpmg.toStdVector(),r2eff.toStdVector(),
                     error.toStdVector(), temp_a, ia, carverrichards, carverrichards);
        calc.fit();
        for (unsigned int i=0; i<a.size(); ++i) a[i].push_back(calc.a[i]);
    }    
}

void Jackknife::compute_cpmg() {
    for (int i=0;i<dataset->n_cpmgVector.size();i++) {
        n_cpmg.resize(0);
        error.resize(0);
        r2eff.resize(0);
        for (int j=0;j<dataset->n_cpmgVector.size();j++)
            if (i!=j) {
                n_cpmg.push_back(dataset->n_cpmgVector[j]);
                r2eff.push_back(dataset->R2effVector[j]);
                error.push_back(dataset->dyVector[j]);
            }
        calculate();
    }
    QVector<double> average (signed (a.size()), 0.);
    double kd_average{};
    double ndata = a[0].size();
    double njack = a[0].size()-1;

    std::vector<double>kd_vector{}, a_average(a.size(), 0.);
    for (unsigned int i=0;i<ndata ;i++) {
        std::vector<double> a_temp{a[0][i],a[1][i],a[2][i],a[3][i],dataset->ligand_cpmg};
        kd_vector.push_back(CPMG_kd(a_temp));
        kd_average+=kd_vector[i];
        for (int j=0; j<signed(a.size()); ++j)
            a_average[j] += a[j][i];

    }
    kd_average/=ndata;
    for (int i=0; i<signed(a_average.size()); ++i) a_average[i] /= ndata;

    double sum{};
    std::vector<double> a_sum(a_average.size(), 0.);
    for (unsigned int i=0;i<ndata;i++){
        sum+=(kd_vector[i]-kd_average)*(kd_vector[i]-kd_average);
        for (int j=0; j<signed(a.size()); ++j)
            a_sum[j] += (a[j][i]-a_average[j])*(a[j][i]-a_average[j]);
    }
    dataset->kd_error=sqrt(njack/ndata * sum);
    for (int i=0; i<signed(a.size()); ++i)
        a_sum[i] = (njack/ndata)*a_sum[i];
}
