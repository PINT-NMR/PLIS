#ifndef CALC_H
#define CALC_H

#include <QVector>
#include <cmath>
#include "data.h"
#include "nr3.h"
#include <complex>

void fit_one_site_dilution(const double x1, const double x2, const std::vector<double> &a, double &y, std::vector<double> &dyda);

double fit_one_site_dilution(const double x1, const double x2, const std::vector<double> &a);

void fit_two_sites_dilution(const double x1, const double x2, const std::vector<double> &a, double &y, std::vector<double> &dyda);

double fit_two_sites_dilution(const double x1, const double x2, const std::vector<double> &a);

void fit_four_sites_dilution(const double Ltot, const double Ptot, const std::vector<double> &a, double &y, std::vector<double> &dyda);

double fit_four_sites_dilution(const double Ltot, const double Ptot, const std::vector<double> &a);

void fit_comp_dilution(const double x1, const double x2, const double x3,const std::vector<double> &a, double &y, std::vector<double> &dyda);

double fit_comp_dilution(const double ligand_conc, const double protein_conc, const double comp_conc, const std::vector<double> &a);

void fit_one_site(const double Ltot, const double Ptot, const std::vector<double> &a, double &y, std::vector<double> &dyda);

double fit_one_site(const double Ltot, const double Ptot, const std::vector<double> &a);

void fit_two_sites(const double ligand_conc, const double protein_conc, const std::vector<double> &a, double &y, std::vector<double> &dyda);

double fit_two_sites(const double ligand_conc, const double protein_conc, const std::vector<double> &a);

void fit_four_sites(const double Ltot, const double Ptot, const std::vector<double> &a, double &y, std::vector<double> &dyda);

double fit_four_sites(const double Ltot, const double Ptot, const std::vector<double> &a);

void fit_comp(const double x1, const double x2, const double x3,const std::vector<double> &a, double &y, std::vector<double> &dyda);

double fit_comp(const double ligand_conc, const double protein_conc, const double comp_conc, const std::vector<double> &a);

double ligand_free_comp(const double start_ligand, const double protein_conc, const double comp_conc, const std::vector<double> &a, const double Lfree);

bool gaussj(MatDoub_IO &a, MatDoub_IO &b);

template <class T>
double zbrent(T &func, const double x1, const double x2, const double tol, const double xx1, const double xx2, const std::vector<double> &aa);
template <class T>
double zbrent2(T &func, const double x1, const double x2, const double tol, const double xx1, const double xx2, const double xx3,const std::vector<double> &aa);

void guess_parameters_one_site(Data* &data, QVector<double> kd_testVector, QString mode);

void guess_parameters_range(Data *&data, QVector<double> &kd, QVector<double> &start, QVector<double> &end, QString mode);

void guess_parameters_range_2_sites(Data *&data, QVector<double> &kd, QVector<double> &kd2, QVector<double> &start, QVector<double> &inter, QVector<double> &end, QString mode);

void guess_parameters_two_sites(Data *&data, QVector<double> kd_testVector, QProgressDialog &progress, QString mode);

void guess_parameters_range_four_sites(Data *&data, QVector<double> &kd1_testVector, QVector<double> &kd2_testVector, QVector<double> &kd3_testVector,
                                       QVector<double> &kd4_testVector, QVector<double> &p_resp_testVector, QVector<double> &pl_resp_testVector,
                                       QVector<double> &pl2_resp_testVector, QVector<double> &pl3_resp_testVector, QVector<double> &pl4_resp_testVector,
                                       QString mode);

void guess_parameters_four_sites(Data *&data, QProgressDialog &progress, QString mode);

void guess_parameters_comp(Data *&data, QVector<double> kd_testVector, QProgressDialog &progress, QString mode);

void guess_parameters_comp_range(Data *&data, QVector<double> &kd, QVector<double> &kd2, QVector<double> &kdc, QVector<double> &start, QVector<double> &end, QString mode);

void guess_parameters_cpmg(Data *&data, QProgressDialog &progress);

void guess_parameters_cpmg_range(Data *&data, QVector<double> &r20, QVector<double> &kex,
                                 QVector<double> &pb, QVector<double> &dw);

//double secord_deriv_dilute(const double ligand, unsigned int idx, const double protein, const std::vector<double> &a,
//                           const double comp=(-1.0));

double RexCR(double R20, double DW, double PA, double KEX, double nuCP);

void carverrichards(double x, VecDoub_I &a, double &y, VecDoub_O &dyda);

double carverrichards(double x, VecDoub_I &a);

inline std::complex<double> acosh_complex(std::complex<double> z);

double CPMG_kd(VecDoub_I &a);

#endif
