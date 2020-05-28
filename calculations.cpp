#include "calculations.h"
#include <iostream>

//One site
void fit_one_site_dilution(const double x1, const double x2, const std::vector<double> &a, double &y, std::vector<double> &dyda){
    double temp{};

    temp=((x1/x2)+(((x2-x1+a[2])/2)-sqrt((((x2-x1+a[2])*(x2-x1+a[2]))/4)+a[2]*x1))/x2);

    y = ((a[0] + (a[1]-a[0])*temp)*x2)/a.back() ;

    dyda[0]= (1 - temp)*x2/a.back();
    dyda[1]= temp*x2/a.back();
    dyda[2]=x2*(a[1]-a[0])*(0.5-((((a[2]+x2-x1)/2)+x1) / sqrt((((a[2]+x2-x1)*(a[2]+x2-x1))/4)+x1*a[2])))/(x2*a.back());

}

double fit_one_site_dilution(const double x1, const double x2, const std::vector<double> &a){
    double temp{};
    double y{};

    temp=((x1/x2)+(((x2-x1+a[2])/2)-sqrt((((x2-x1+a[2])*(x2-x1+a[2]))/4)+a[2]*x1)) /x2);

    y = (a[0] + (a[1]-a[0])*temp)*x2/a.back();
    return y;
}

//Two-site
void fit_two_sites_dilution(const double ligand_conc, const double protein_conc, const std::vector<double> &a, double &y, std::vector<double> &dyda){
    double P{}, root{}, start{a[a.size()-2]}; //a[5] is start conc of the datasets/conc_vector, since we need it to make zbrent

    try {
        root=zbrent(ligand_free,start,ligand_conc,1e-12, ligand_conc,protein_conc,a);
    } catch (int) {
        //std::cout << "heeej" << std::endl;
    }
    P=protein_free(root,protein_conc,a);
    double PL{}, PL2{};
    PL=P*root/a[3];
    PL2=PL*root/a[4];
    y = ((a[0]*P + a[1]*PL + a[2] * PL2)/protein_conc)*protein_conc/a.back();

    for (unsigned int i=0;i<dyda.size();i++)
        dyda[i]=secord_deriv_dilute(ligand_conc, i, protein_conc, a);
}

double fit_two_sites_dilution(const double ligand_conc, const double protein_conc, const std::vector<double> &a){
    double P{}, root{}, y{}, start{a[a.size()-2]}; //a[5] is start conc of the datasets/conc_vector, since we need it to make zbrent
    try {
        root=zbrent(ligand_free,start,ligand_conc,1e-12, ligand_conc,protein_conc,a);
    } catch (int) {

    }
    P=protein_free(root,protein_conc,a);
    double PL{}, PL2{};
    PL=P*root/a[3];
    PL2=PL*root/a[4];
    y = ((a[0]*P + a[1]*PL + a[2] * PL2)/protein_conc)*protein_conc/a.back();
    return y;
}

double ligand_free(const double start_ligand, const double protein_conc, const std::vector<double> &a, const double Lfree){
    return (Lfree*Lfree*Lfree/(a[3]*a[4]))+(Lfree*Lfree)*((1/a[3])+(2*protein_conc/(a[3]*a[4]))-1/(a[3]*a[4]))+
            Lfree*(1+(protein_conc/a[3])-(start_ligand/a[3]))-start_ligand;
}

double protein_free(const double L, const double x2, const std::vector<double> &a){
    return x2/(1+(L/a[3])+(L*L)/(a[3]*a[4]));
}


//Competitive Binding
void fit_comp(const double ligand_conc, const double protein_conc, const double comp_conc, const std::vector<double> &a, double &y, std::vector<double> &dyda){
    double I{},IL{}, root{}, start{a[a.size()-3]}; //a[5] is start conc of the datasets/conc_vector, since we need it to make zbrent
    try {
        root=zbrent2(ligand_free_comp,start,ligand_conc,1e-12, ligand_conc,protein_conc,comp_conc,a);
    } catch (int) {

    }
    I=(a[4]*comp_conc)/(a[4]+root);
    IL=comp_conc-I;
    y = ((a[0]*I + a[1]*IL)/comp_conc)*comp_conc/a.back();

    for (unsigned int i=0;i<dyda.size();i++)
        dyda[i]=secord_deriv_dilute(ligand_conc, i, protein_conc, a, comp_conc);
}

double fit_comp(const double ligand_conc, const double protein_conc, const double comp_conc, const std::vector<double> &a){
    double I{},IL{}, y{}, root{}, start{a[a.size()-3]}; //a[5] is start conc of the datasets/conc_vector, since we need it to make zbrent
    try {
        root=zbrent2(ligand_free_comp,start,ligand_conc,1e-12, ligand_conc,protein_conc,comp_conc,a);
    } catch (int) {

    }
    I=(a[4]*comp_conc)/(a[4]+root);
    IL=comp_conc-I;
    y = ((a[0]*I + a[1]*IL)/comp_conc)*comp_conc/a.back();
    return y;
}

double ligand_free_comp(const double start_ligand, const double protein_conc, const double comp_conc, const std::vector<double> &a, const double Lfree){
    return Lfree+comp_conc*Lfree/(a[4]+Lfree)+(protein_conc/(1+Lfree/a[2]+2*(Lfree*Lfree)/(a[2]*a[3])))*(Lfree/a[2]+2*(Lfree*Lfree)/(a[2]*a[3]))-start_ligand;
}

bool gaussj(MatDoub_IO &a, MatDoub_IO &b)
{
    Int i,icol,irow,j,k,l,ll,n=a.nrows(),m=b.ncols();
    Doub big,dum,pivinv;
    VecInt indxc(n),indxr(n),ipiv(n);
    for (j=0;j<n;j++) ipiv[j]=0;
    for (i=0;i<n;i++) {
        big=0.0;
        for (j=0;j<n;j++)
            if (ipiv[j] != 1)
                for (k=0;k<n;k++) {
                    if (ipiv[k] == 0) {
                        if (abs(a[j][k]) >= big) {
                            big=abs(a[j][k]);
                            irow=j;
                            icol=k;
                        }
                    }
                }
        ++(ipiv[icol]);
        if (irow != icol) {
            for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
            for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l]);
        }
        indxr[i]=irow;
        indxc[i]=icol;
        if (a[icol][icol] == 0.0) return false;//throw("gaussj: Singular Matrix");
        pivinv=1.0/a[icol][icol];
        a[icol][icol]=1.0;
        for (l=0;l<n;l++) a[icol][l] *= pivinv;
        for (l=0;l<m;l++) b[icol][l] *= pivinv;
        for (ll=0;ll<n;ll++)
            if (ll != icol) {
                dum=a[ll][icol];
                a[ll][icol]=0.0;
                for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
                for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
            }
    }
    for (l=n-1;l>=0;l--) {
        if (indxr[l] != indxc[l])
            for (k=0;k<n;k++)
                SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    }
    return true;
}

template <class T>
Doub zbrent(T &func, const Doub x1, const Doub x2, const Doub tol, const double xx1, const double xx2, const std::vector<double> &aa)
{
    const Int ITMAX=100;
    const Doub EPS=numeric_limits<Doub>::epsilon();
    Doub a=x1,b=x2,c=x2,d,e,fa=func(xx1,xx2, aa, a),fb=func(xx1,xx2, aa, b),fc,p,q,r,s,tol1,xm;
    //std::cout << "fa = " << fa << ", fb = " << fb <<  ", x1 = "<<x1 << ", x2 = " << x2 << std::endl;
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        throw("Root must be bracketed in zbrent");
    fc=fb;
    for (Int iter=0;iter<ITMAX;iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        if (abs(fc) < abs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*EPS*abs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (abs(xm) <= tol1 || fb == 0.0) return b;
        if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q;
            p=abs(p);
            Doub min1=3.0*xm*q-abs(tol1*q);
            Doub min2=abs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {
                e=d;
                d=p/q;
            } else {
                d=xm;
                e=d;
            }
        } else {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if (abs(d) > tol1)
            b += d;
        else
            b += SIGN(tol1,xm);
            fb=func(xx1,xx2, aa, b);
    }
    throw("Maximum number of iterations exceeded in zbrent");
}

template <class T>
Doub zbrent2(T &func, const Doub x1, const Doub x2, const Doub tol, const double xx1, const double xx2, const double xx3,const std::vector<double> &aa)
{
    const Int ITMAX=100;
    const Doub EPS=numeric_limits<Doub>::epsilon();
    Doub a=x1,b=x2,c=x2,d,e,fa=func(xx1,xx2,xx3, aa, a),fb=func(xx1,xx2,xx3, aa, b),fc,p,q,r,s,tol1,xm;
    //std::cout << "fa = " << fa << ", fb = " << fb <<  ", x1 = "<<x1 << ", x2 = " << x2 << std::endl;
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        throw("Root must be bracketed in zbrent");
    fc=fb;
    for (Int iter=0;iter<ITMAX;iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        if (abs(fc) < abs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*EPS*abs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (abs(xm) <= tol1 || fb == 0.0) return b;
        if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q;
            p=abs(p);
            Doub min1=3.0*xm*q-abs(tol1*q);
            Doub min2=abs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {
                e=d;
                d=p/q;
            } else {
                d=xm;
                e=d;
            }
        } else {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if (abs(d) > tol1)
            b += d;
        else
            b += SIGN(tol1,xm);
            fb=func(xx1,xx2, xx3,aa, b);
    }
    throw("Maximum number of iterations exceeded in zbrent");
}

//All guess parameter functions is when the software estimates a starting guess
//All guess parameters_range is when the user has guessed within a range.

void guess_parameters_one_site(Data *&data, QVector<double> kd_testVector){
    double y_calc{};
    double chi2{-1.0};
    double chi2_temp{};
    int best_index_kd{}, best_index_responce{};
    QVector<double> test_last_responce{};
    double range{data->responceVector.last()/10};
    for (int i=0;i<11;i++)
        test_last_responce.push_back(data->responceVector.last()+i*range);

    for (int i=0; i<kd_testVector.size();i++){
        std::vector<double> a_temp{data->a};
        a_temp[2]=kd_testVector[i];
        for (int h=0;h<test_last_responce.size();h++){
            a_temp[1]=test_last_responce[h];
            for (int i=0; i<data->responceVector.size();i++){
                y_calc=fit_one_site_dilution(data->concVector[i],data->protein_conc_vector[i], a_temp);
                chi2_temp += ((y_calc - data->responceVector[i]) *
                              (y_calc - data->responceVector[i]))
                        / (data->error_vector[i]*data->error_vector[i]);
            }

            if(!isnan(chi2_temp)){
                if (chi2<0){
                    chi2=chi2_temp;
                    best_index_kd=i;
                    best_index_responce=h;
                }
                else if(chi2>chi2_temp){
                    chi2=chi2_temp;
                    best_index_kd=i;
                    best_index_responce=h;
                }
            }
            chi2_temp=0;
        }
    }
    data->a[2]=kd_testVector[best_index_kd];
    data->a[1]=test_last_responce[best_index_responce];
}

void guess_parameters_range(Data *&data, QVector<double> &kd, QVector<double> &start, QVector<double> &end){
    double y_calc{};
    double chi2{};
    double chi2_temp{};
    int best_index_kd{}, best_index_start_responce{}, best_index_end_responce{};
    std::vector<double> a_temp{double(),double(),double()};
    for (int i = 0; i < kd.size();i++) {
        a_temp[2]=kd[i];
        for (int j=0; j < start.size();j++){
            a_temp[0]=start[j];
            for (int l=0;l<end.size();l++){
                a_temp[1]=end[l];
                for (int k=0; k<data->responceVector.size();k++){
                    y_calc=fit_one_site_dilution(data->concVector[k],data->protein_conc_vector[k], a_temp);
                    chi2_temp += ((y_calc - data->responceVector[k]) *
                                  (y_calc - data->responceVector[k]))
                            / (data->error_vector[k]*data->error_vector[k]);
                }
                if (chi2==0.0 && i==0 && j==0 && l==0){
                    chi2=chi2_temp;
                    best_index_end_responce=l;
                    best_index_start_responce=j;
                    best_index_kd=i;
                }
                else if(chi2>chi2_temp){
                    chi2=chi2_temp;
                    best_index_end_responce=l;
                    best_index_start_responce=j;
                    best_index_kd=i;
                }
                chi2_temp=0;
            }
        }
    }
    data->a[0]=start[best_index_start_responce];
    data->a[1]=end[best_index_end_responce];
    data->a[2]=kd[best_index_kd];
}

void guess_parameters_range_2_sites(Data *&data, QVector<double> &kd, QVector<double> &kd2, QVector<double> &start, QVector<double> &inter, QVector<double> &end){
    double y_calc{};
    double chi2{-1.0};
    double chi2_temp{};
    int best_index_kd{}, best_index_start_responce{}, best_index_end_responce{}, best_index_inter_responce{}, best_index_kd2{};
    std::vector<double> a_temp{data->a};
    for (int i = 0; i < kd.size();i++) {
        a_temp[3]=kd[i];
        for (int j=0; j < start.size();j++){
            a_temp[0]=start[j];
            for (int l=0; l < end.size();l++){
                a_temp[2]=end[l];
                for (int m=0;m<inter.size();m++){
                    a_temp[1]=inter[m];
                    for (int n=0;n<kd2.size();n++){
                        a_temp[4]=kd2[n];
                        for (int k=0; k<data->responceVector.size();k++){
                            y_calc=fit_two_sites_dilution(data->concVector[k],data->protein_conc_vector[k], a_temp);
                            chi2_temp += ((y_calc - data->responceVector[k]) *
                                          (y_calc - data->responceVector[k]))
                                    / (data->error_vector[k]*data->error_vector[k]);
                        }
                        if(!isnan(chi2_temp)){
                            if (chi2<0){
                                chi2=chi2_temp;
                                best_index_end_responce=l;
                                best_index_start_responce=j;
                                best_index_kd=i;
                                best_index_inter_responce=m;
                                best_index_kd2=n;
                            }
                            else if(chi2>chi2_temp){
                                chi2=chi2_temp;
                                best_index_end_responce=l;
                                best_index_start_responce=j;
                                best_index_kd=i;
                                best_index_inter_responce=m;
                                best_index_kd2=n;
                            }
                        }
                        chi2_temp=0;
                    }
                }
            }
        }
    }
    data->a[0]=start[best_index_start_responce];
    data->a[1]=inter[best_index_inter_responce];
    data->a[2]=end[best_index_end_responce];
    data->a[3]=kd[best_index_kd];
    data->a[4]=kd2[best_index_kd2];
}

void guess_parameters_two_sites(Data *&data, QVector<double> kd_testVector, QProgressDialog &progress){
    double y_calc{};
    double chi2_temp{};
    double chi2{-1.0};

    int best_index_kd1{}, best_index_kd2{};
    QVector<double> test_last_responce{}, test_intermediate_responce{};

    double range{data->responceVector.last()/10};
    for (int i=0;i<11;i++)
        test_last_responce.push_back(data->responceVector.last()+i*range);
    for (int i=0;i<data->responceVector.size();i++)
        test_intermediate_responce.push_back(data->responceVector[i]);
    for (int i=1;i<11;i++)
        test_intermediate_responce.push_back(test_last_responce[i]);
    int best_index_a1{}, best_index_a2{};

    std::vector<double> a_temp{data->a};
    for (int g = 0; g < test_last_responce.size();g++) {
        progress.setValue(g+1);
        a_temp[2]=test_last_responce[g];
        for (int h=0; h < test_intermediate_responce.size();h++){
            a_temp[1]=test_intermediate_responce[h];
            for (int i=0; i<kd_testVector.size();i++){
                a_temp[3]=kd_testVector[i];
                for (int j=0;j<kd_testVector.size();j++){
                    a_temp[4]=kd_testVector[j];
                    for (int k=0; k<data->responceVector.size();k++){
                        y_calc=fit_two_sites_dilution(data->concVector[k],data->protein_conc_vector[k], a_temp);
                        chi2_temp += ((y_calc - data->responceVector[k]) *
                                      (y_calc - data->responceVector[k]))
                                / (data->error_vector[k]*data->error_vector[k]);
                    }
                    if(!isnan(chi2_temp)){
                        if (chi2<0){
                            chi2=chi2_temp;
                            best_index_kd1=i;
                            best_index_kd2=j;
                            best_index_a1=h;
                            best_index_a2=g;
                        }
                        else if(chi2>chi2_temp){
                            chi2=chi2_temp;
                            best_index_kd1=i;
                            best_index_kd2=j;
                            best_index_a1=h;
                            best_index_a2=g;
                        }
                    }
                    chi2_temp=0;
                }
            }
        }
    }
    progress.setValue(test_last_responce.size()+1);
    data->a[3]=kd_testVector[best_index_kd1];
    data->a[4]=kd_testVector[best_index_kd2];
    data->a[1]=test_intermediate_responce[best_index_a1];
    data->a[2]=test_last_responce[best_index_a2];
}

void guess_parameters_comp(Data *&data, QVector<double> kd_testVector, QProgressDialog &progress){
    double y_calc{};
    double chi2_temp{};
    double chi2{(-1.0)};
    int best_index_kd1{}, best_index_kd2{}, best_index_kdc{},best_index_a1{};

    QVector<double> test_last_responce{};
    double range{data->responceVector.last()/10};
    if(data->responceVector.back()>data->responceVector[0]){
        for (int i=0;i<11;i++)
            test_last_responce.push_back(data->responceVector.last()+i*range);
    }
    else
        for (int i=0;i<11;i++)
            test_last_responce.push_back(data->responceVector.last()-i*range);

    std::vector<double> a_temp{data->a};
    for (int g = 0; g < test_last_responce.size();g++) {
        a_temp[1]=test_last_responce[g];
        progress.setValue(g+1);
        for (int i=0; i<kd_testVector.size();i++){
            a_temp[2]=kd_testVector[i];
            for (int j=i;j<kd_testVector.size();j++){
                a_temp[3]=kd_testVector[j];
                for (int c=0;c<kd_testVector.size();c++){
                    a_temp[4]=kd_testVector[c];
                    for (int k=0; k<data->responceVector.size();k++){
                        y_calc=fit_comp(data->concVector[k],data->protein_conc_vector[k],data->comp_vector[k], a_temp);
                        chi2_temp += ((y_calc - data->responceVector[k]) *
                                      (y_calc - data->responceVector[k]))
                                / (data->error_vector[k]*data->error_vector[k]);
                    }
                    if(!isnan(chi2_temp)){
                        if (chi2<0){
                            chi2=chi2_temp;
                            best_index_kd1=i;
                            best_index_kd2=j;
                            best_index_a1=g;
                            best_index_kdc=c;
                        }
                        else if(chi2>chi2_temp){
                            chi2=chi2_temp;
                            best_index_kd1=i;
                            best_index_kd2=j;
                            best_index_a1=g;
                            best_index_kdc=c;
                        }
                    }
                    chi2_temp=0;
                }
            }
        }
    }
    progress.setValue(test_last_responce.size()+1);
    data->a[2]=kd_testVector[best_index_kd1];
    data->a[3]=kd_testVector[best_index_kd2];
    data->a[4]=kd_testVector[best_index_kdc];
    data->a[1]=test_last_responce[best_index_a1];
}

void guess_parameters_comp_range(Data *&data, QVector<double> &kd, QVector<double> &kd2, QVector<double> &kdc, QVector<double> &start, QVector<double> &end){
    double y_calc{};
    double chi2{-1.0};
    double chi2_temp{};
    int best_index_kd1{}, best_index_start_responce{}, best_index_end_responce{}, best_index_kdc{}, best_index_kd2{};
    std::vector<double> a_temp{data->a};
    for (int i = 0; i < kd.size();i++) {
        a_temp[2]=kd[i];
        for (int j=0; j < start.size();j++){
            a_temp[0]=start[j];
            for (int l=0; l < end.size();l++){
                a_temp[1]=end[l];
                for (int m=0;m<kdc.size();m++){
                    a_temp[4]=kdc[m];
                    for (int n=0;n<kd2.size();n++){
                        a_temp[3]=kd2[n];
                        for (int k=0; k<data->responceVector.size();k++){
                            y_calc=fit_comp(data->concVector[k],data->protein_conc_vector[k],data->comp_vector[k], a_temp);
                            chi2_temp += ((y_calc - data->responceVector[k]) *
                                          (y_calc - data->responceVector[k]))
                                    / (data->error_vector[k]*data->error_vector[k]);
                        }
                        if(!isnan(chi2_temp)){
                            if (chi2<0){
                                chi2=chi2_temp;
                                best_index_end_responce=l;
                                best_index_start_responce=j;
                                best_index_kd1=i;
                                best_index_kdc=m;
                                best_index_kd2=n;
                            }
                            else if(chi2>chi2_temp){
                                chi2=chi2_temp;
                                best_index_end_responce=l;
                                best_index_start_responce=j;
                                best_index_kd1=i;
                                best_index_kdc=m;
                                best_index_kd2=n;
                            }
                        }
                        chi2_temp=0;
                    }
                }
            }
        }
    }
    data->a[0]=start[best_index_start_responce];
    data->a[1]=end[best_index_end_responce];
    data->a[2]=kd[best_index_kd1];
    data->a[3]=kd2[best_index_kd2];
    data->a[4]=kdc[best_index_kdc];
}

void guess_parameters_cpmg(Data *&data, QProgressDialog &progress){
    double y_calc{};
    double chi2_temp{};
    double chi2{};

    int best_index_r20{}, best_index_pb{}, best_index_kex{}, best_index_dw{};
    QVector<double> r20_test{}, dw_test{},kex_test{}, pb_test{0.005, 0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 0.30, 0.50};

    for (int i=0;i<11;i++)
        r20_test.push_back(2.5*i+5.0);
    for (int i=0;i<11;i++)
        dw_test.push_back(450*i+500.0);
    for (int i=0;i<11;i++)
        kex_test.push_back(190*i+100.0);

    std::vector<double> a_temp{data->a};
    for (int g = 0; g < r20_test.size();g++) {
        progress.setValue(g+1);
        a_temp[0]=r20_test[g];
        for (int h=0; h < kex_test.size();h++){
            a_temp[1]=kex_test[h];
            for (int i=0; i<pb_test.size();i++){
                a_temp[2]=pb_test[i];
                for (int j=i;j<dw_test.size();j++){
                    a_temp[3]=dw_test[j];
                    for (int k=0; k<data->n_cpmgVector.size();k++){
                        y_calc=carverrichards(data->n_cpmgVector[k], a_temp);
                        chi2_temp += ((y_calc - data->R2effVector[k]) *
                                      (y_calc - data->R2effVector[k]))
                                / (data->dyVector[k]*data->dyVector[k]);
                    }
                    if(!isnan(chi2_temp)){
                        if (chi2<0){
                            chi2=chi2_temp;
                            best_index_dw=j;
                            best_index_r20=g;
                            best_index_pb=i;
                            best_index_kex=h;
                        }
                        else if(chi2>chi2_temp){
                            chi2=chi2_temp;
                            best_index_dw=j;
                            best_index_r20=g;
                            best_index_pb=i;
                            best_index_kex=h;
                        }
                    }
                    chi2_temp=0;
                }
            }
        }
    }
    progress.setValue(r20_test.size()+1);
    data->a[0]=r20_test[best_index_r20];
    data->a[1]=kex_test[best_index_kex];
    data->a[2]=pb_test[best_index_pb];
    data->a[3]=dw_test[best_index_dw];
}


void guess_parameters_cpmg_range(Data *&data, QVector<double> &r20, QVector<double> &kex, QVector<double> &pb, QVector<double> &dw){
    double y_calc{};
    double chi2{-1.0};
    double chi2_temp{};
    int best_index_r20{}, best_index_pb{}, best_index_kex{}, best_index_dw{};
    std::vector<double> a_temp{data->a};
    for (int i = 0; i < r20.size();i++) {
        a_temp[0]=r20[i];
        for (int j=0; j < dw.size();j++){
            a_temp[3]=dw[j];
            for (int m=0;m<pb.size();m++){
                a_temp[2]=pb[m];
                for (int n=0;n<kex.size();n++){
                    a_temp[1]=kex[n];
                    for (int k=0; k<data->n_cpmgVector.size();k++){
                        y_calc=carverrichards(data->n_cpmgVector[k], a_temp);
                        chi2_temp += ((y_calc - data->R2effVector[k]) *
                                      (y_calc - data->R2effVector[k]))
                                / (data->dyVector[k]*data->dyVector[k]);
                    }
                    if(!isnan(chi2_temp)){
                        if (chi2<0){
                            chi2=chi2_temp;
                            best_index_dw=j;
                            best_index_r20=i;
                            best_index_pb=m;
                            best_index_kex=n;
                        }
                        else if(chi2>chi2_temp){
                            chi2=chi2_temp;
                            best_index_dw=j;
                            best_index_r20=i;
                            best_index_pb=m;
                            best_index_kex=n;
                        }
                    }
                    chi2_temp=0;
                }
            }
        }
    }
    data->a[0]=r20[best_index_r20];
    data->a[1]=kex[best_index_kex];
    data->a[2]=pb[best_index_pb];
    data->a[3]=dw[best_index_dw];
}

double secord_deriv_dilute(const double ligand, unsigned int idx, const double protein, const std::vector<double> &a, const double comp)
{
// Error bound ~eps^2/3
    std::vector<double> a_temp{a};
    const double eps = std::numeric_limits<double>::epsilon();
    double aidx = a_temp[idx];
    double h = pow(3*eps, 0.33333333333333);
    volatile double temp = a_temp[idx] + h;
    h = temp - a_temp[idx];
    if (h == 0.0) { h = nextafter(a_temp[idx], std::numeric_limits<double>::max()) - a_temp[idx]; }
    a_temp[idx] = aidx + h;
    double yh{},ymh{};
    if(comp==(-1.0)){
        yh = fit_two_sites_dilution(ligand,protein, a_temp);
        a_temp[idx] = aidx - h;
        ymh = fit_two_sites_dilution(ligand,protein, a_temp);
    }
    else{ //befor, fit_two_sites_diluted, not comp, borde vara comp? med bruset på 0,002 så blire kaos, så kan man inte göra anpassningen
        yh = fit_comp(ligand,protein,comp, a_temp);
        a_temp[idx] = aidx - h;
        ymh = fit_comp(ligand,protein,comp, a_temp);
    }
    double diff = yh - ymh;
    a_temp[idx] = aidx;
    return diff/(2*h);
}

Doub RexCR(Doub R20, Doub DW, Doub PB, Doub KEX, Doub nuCP)
{
// R20 = (R2A + R2B)/2
// Note that order of arguments does not correspond to order of a-vector!!
// Also note that 2*tauCP is separation of 180-pulses

    complex<Doub> psi, zeta, np, nm, Dp, Dm;
    Doub tmp;
    Doub PA = 1.-PB;
    Doub tauCP = 0.25/nuCP;

    psi = SQR(-PA*KEX+PB*KEX)-SQR(DW)+4.*PA*PB*SQR(KEX);
    zeta = 2.*DW*(-PA*KEX+PB*KEX);
    np = sqrt(2.)*tauCP*sqrt(sqrt(SQR(psi) + SQR(zeta)) + psi);
    nm = sqrt(2.)*tauCP*sqrt(sqrt(SQR(psi) + SQR(zeta)) - psi);
    Dp = 0.5*((psi + 2.*SQR(DW))/sqrt(SQR(psi) + SQR(zeta)) + 1.);
    Dm = 0.5*((psi + 2.*SQR(DW))/sqrt(SQR(psi) + SQR(zeta)) - 1.);
    tmp = fabs(real(acosh_complex(Dp*cosh(np)-Dm*cos(nm))));
    return 0.5*(2.*R20+KEX-1./(2.*tauCP)*tmp);
}

void carverrichards(Doub x, VecDoub_I &a, Doub &y, VecDoub_O &dyda)
{
    const Doub STEP = 1.0e-6;
    Doub R20 = a[0];
    Doub KEX = a[1];
    Doub PB = a[2];
    Doub DW = a[3];

    y = RexCR(R20, DW, PB, KEX, x);

    dyda[0] = 1.0;
    dyda[1] = (RexCR(R20, DW, PB, (1.+STEP)*KEX, x) - y ) / (STEP*KEX);
    dyda[2] = (RexCR(R20, DW, (1.+STEP)*PB, KEX, x) - y ) / (STEP*PB);
    dyda[3] = (RexCR(R20, (1.+STEP)*DW, PB, KEX, x) - y ) / (STEP*DW);
}

double carverrichards(Doub x, VecDoub_I &a)
{
    Doub R20 = a[0];
    Doub KEX = a[1];
    Doub PB = a[2];
    Doub DW = a[3];

    return RexCR(R20, DW, PB, KEX, x);
}

inline complex<Doub> acosh_complex(complex<Doub> z) {
    return log(z + sqrt(z-1.)*sqrt(z+1.));
}

double CPMG_kd(VecDoub_I &a){
    Doub PB = a[2];
    Doub Kb = a[1]-PB*a[1];
    Doub Kf_L = a[1]-Kb;
    return (Kb * a.back())/ Kf_L;
}
