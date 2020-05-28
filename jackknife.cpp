#include "jackknife.h"

Jackknife::Jackknife(Data *&data)
    :dataset(data)
{

}

void Jackknife::compute(){
    for (int i=0;i<dataset->responceVector.size();i++){
        responce.clear();
        conc.clear();
        protein.clear();
        error.clear();
        volume.clear();
        comp.clear();
        if (dataset->num_bind_site==1 || dataset->num_bind_site==2)
            for (int j=0;j<dataset->responceVector.size();j++){
                if (i!=j){
                    responce.push_back(dataset->responceVector[j]);
                    conc.push_back(dataset->concVector[j]);
                    protein.push_back(dataset->protein_conc_vector[j]);
                    error.push_back(dataset->error_vector[j]);
                    volume.push_back(dataset->volumeVector[j]);
                }
            }
        else if(dataset->num_bind_site==3)
            for (int j=0;j<dataset->responceVector.size();j++){
                if (i!=j){
                    responce.push_back(dataset->responceVector[j]);
                    conc.push_back(dataset->concVector[j]);
                    protein.push_back(dataset->protein_conc_vector[j]);
                    error.push_back(dataset->error_vector[j]);
                    volume.push_back(dataset->volumeVector[j]);
                    comp.push_back(dataset->comp_vector[j]);
                }
            }
        calculate();
    }
    QVector<double> average (signed (a.size()),0.);
    double ndata;
    double jack;
    jack=a[0].size()-1;
    ndata=a[0].size();
    int loop{};
    if(dataset->num_bind_site==1)
        loop=3;
    else if(dataset->num_bind_site==2 || dataset->num_bind_site==3)
        loop=5;

    for(unsigned int j=0;j<unsigned (loop);j++){
        for (unsigned int k=0;k<ndata;k++)
            average[signed (j)]+=a[j][k];
    }
    if(dataset->num_bind_site==1){
        average[2]/=ndata;
        double sum{};
        for (unsigned int i=0;i<ndata;i++){
            sum+=(a[2][i]-average[2])*(a[2][i]-average[2]);
        }
        dataset->kd_error=sqrt(jack/ndata * sum);
    }
    else if(dataset->num_bind_site==2){
        average[3]/=ndata;
        average[4]/=ndata;
        double sum{},sum2{};
        for (unsigned int i=0;i<ndata;i++){
            sum+=(a[3][i]-average[3])*(a[3][i]-average[3]);
            sum2+=(a[4][i]-average[4])*(a[4][i]-average[4]);
        }
        dataset->kd1_error=sqrt(jack/ndata * sum);
        dataset->kd2_error=sqrt(jack/ndata * sum2);
    }
    else if(dataset->num_bind_site==3){
        average[2]/=ndata;
        average[3]/=ndata;
        average[4]/=ndata;
        double sum{},sum2{}, sum3{};
        for (unsigned int i=0;i<ndata;i++){
            sum+=(a[2][i]-average[2])*(a[2][i]-average[2]);
            sum2+=(a[3][i]-average[3])*(a[3][i]-average[3]);
            sum3+=(a[4][i]-average[4])*(a[4][i]-average[4]);
        }
        dataset->kd1_error=sqrt(jack/ndata * sum);
        dataset->kd2_error=sqrt(jack/ndata * sum2);
        dataset->kdc_error=sqrt(jack/ndata * sum3);
    }
}

void Jackknife::calculate(){
    std::vector<double> temp_a{dataset->a};
    if (dataset->num_bind_site==1){
        Fitmrq calc(conc.toStdVector(),protein.toStdVector(),responce.toStdVector(),
                    volume.toStdVector(),error.toStdVector(),temp_a,fit_one_site_dilution,fit_one_site_dilution);
        calc.ia[4]=false;
        calc.ia[3]=false;
        calc.fit();
        a[0].push_back(calc.a[0]);
        a[1].push_back(calc.a[1]);
        a[2].push_back(calc.a[2]);
    }
    else if(dataset->num_bind_site==2){
        Fitmrq calc(conc.toStdVector(),protein.toStdVector(),responce.toStdVector(),
                    volume.toStdVector(),error.toStdVector(),temp_a,fit_two_sites_dilution, fit_two_sites_dilution);
        calc.ia[5]=false;
        calc.ia[6]=false;
        calc.fit();
        a[0].push_back(calc.a[0]);
        a[1].push_back(calc.a[1]);
        a[2].push_back(calc.a[2]);
        a[3].push_back(calc.a[3]);
        a[4].push_back(calc.a[4]);
    }
    else if(dataset->num_bind_site==3){
        Fitmrq2 calc(conc.toStdVector(),protein.toStdVector(),comp.toStdVector(),responce.toStdVector(),
                    volume.toStdVector(),error.toStdVector(),temp_a,fit_comp,fit_comp);
        calc.ia[5]=false;
        calc.ia[6]=false;
        calc.ia[7]=false;
        calc.fit();
        a[0].push_back(calc.a[0]);
        a[1].push_back(calc.a[1]);
        a[2].push_back(calc.a[2]);
        a[3].push_back(calc.a[3]);
        a[4].push_back(calc.a[4]);
    }
    else if(dataset->num_bind_site==4){
        Fitmrq3 calc(n_cpmg.toStdVector(),r2eff.toStdVector(),
                     error.toStdVector(), temp_a, carverrichards,carverrichards);
        calc.ia[4]=false;
        calc.fit();
        a[0].push_back(calc.a[0]);
        a[1].push_back(calc.a[1]);
        a[2].push_back(calc.a[2]);
        a[3].push_back(calc.a[3]);
        a[4].push_back(calc.a[4]);
    }
}

void Jackknife::compute_cpmg(){
    for (int i=0;i<dataset->n_cpmgVector.size();i++){
        n_cpmg.clear();
        error.clear();
        r2eff.clear();
        for (int j=0;j<dataset->n_cpmgVector.size();j++){
            if (i!=j){
                n_cpmg.push_back(dataset->n_cpmgVector[j]);
                r2eff.push_back(dataset->R2effVector[j]);
                error.push_back(dataset->dyVector[j]);
            }
        }
        calculate();
    }
    QVector<double> average (signed (a.size()),0.);
    double kd_average{};
    double ndata;
    double jack;
    jack=a[0].size()-1;
    ndata=a[0].size();

    std::vector<double>kd_vector{};
    for(unsigned int i=0;i<ndata ;i++){
        std::vector<double> a_temp{a[0][i],a[1][i],a[2][i],a[3][i],dataset->ligand_cpmg};
        kd_vector.push_back(CPMG_kd(a_temp));
        kd_average+=kd_vector[i];
    }
    kd_average/=ndata;

    double sum{};
    for (unsigned int i=0;i<ndata;i++){
        sum+=(kd_vector[i]-kd_average)*(kd_vector[i]-kd_average);
    }
    dataset->kd_error=sqrt(jack/ndata * sum);
}
