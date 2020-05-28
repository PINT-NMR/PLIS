#ifndef CALC_H
#define CALC_H

#include <QVector>
#include <cmath>
#include "data.h"
#include "nr3.h"
#include <complex>

double ligand_free(const double x1, const double x2, const std::vector<double> &a, const double L);

double protein_free(const double L, const double x2, const std::vector<double> &a);

void fit_one_site_dilution(const double x1, const double x2, const std::vector<double> &a, double &y, std::vector<double> &dyda);

double fit_one_site_dilution(const double x1, const double x2, const std::vector<double> &a);

void fit_two_sites_dilution(const double x1, const double x2, const std::vector<double> &a, double &y, std::vector<double> &dyda);

double fit_two_sites_dilution(const double x1, const double x2, const std::vector<double> &a);

void fit_comp(const double x1, const double x2, const double x3,const std::vector<double> &a, double &y, std::vector<double> &dyda);

double fit_comp(const double ligand_conc, const double protein_conc, const double comp_conc, const std::vector<double> &a);

double ligand_free_comp(const double start_ligand, const double protein_conc, const double comp_conc, const std::vector<double> &a, const double Lfree);

bool gaussj(MatDoub_IO &a, MatDoub_IO &b);

template <class T>
double zbrent(T &func, const double x1, const double x2, const double tol, const double xx1, const double xx2, const std::vector<double> &aa);
template <class T>
double zbrent2(T &func, const Doub x1, const Doub x2, const Doub tol, const double xx1, const double xx2, const double xx3,const std::vector<double> &aa);

void guess_parameters_one_site(Data* &data, QVector<double> kd_testVector);

void guess_parameters_range(Data *&data, QVector<double> &kd, QVector<double> &start, QVector<double> &end);

void guess_parameters_range_2_sites(Data *&data, QVector<double> &kd, QVector<double> &kd2, QVector<double> &start, QVector<double> &inter, QVector<double> &end);

void guess_parameters_two_sites(Data *&data, QVector<double> kd_testVector, QProgressDialog &progress);

void guess_parameters_comp(Data *&data, QVector<double> kd_testVector, QProgressDialog &progress);

void guess_parameters_comp_range(Data *&data, QVector<double> &kd, QVector<double> &kd2, QVector<double> &kdc, QVector<double> &start, QVector<double> &end);

void guess_parameters_cpmg(Data *&data, QProgressDialog &progress);

void guess_parameters_cpmg_range(Data *&data, QVector<double> &r20, QVector<double> &kex,
                                 QVector<double> &pb, QVector<double> &dw);

double secord_deriv_dilute(const double ligand, unsigned int idx, const double protein, const std::vector<double> &a,
                           const double comp=(-1.0));

Doub RexCR(Doub R20, Doub DW, Doub PA, Doub KEX, Doub nuCP);

void carverrichards(Doub x, VecDoub_I &a, Doub &y, VecDoub_O &dyda);

double carverrichards(Doub x, VecDoub_I &a);

inline complex<Doub> acosh_complex(complex<Doub> z);

double CPMG_kd(VecDoub_I &a);

#endif
