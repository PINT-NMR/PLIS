#include "calculations.h"
#include <iostream>

// 200612
// - took absolute values of all fitting parameters that correspond to Kd:s
// - N.B. This should never be done when calculating the derivatives!!

// ********************************************************************************************
// Dilute means that total concentration must be taken into account for calculating the signal
// *********************************************************************************************

// One site
// -------
// a0-a1: P and PL signals
// a2: Kd
// a3: [L]start
// a4: [P]start

void fit_one_site_dilution(const double Ltot, const double Ptot, const std::vector<double> &a, double &y, std::vector<double> &dyda)
{
    double temp=((Ltot/Ptot)+(((Ptot-Ltot+fabs(a[2]))/2) -
                 sqrt((((Ptot-Ltot+fabs(a[2]))*(Ptot-Ltot+fabs(a[2])))/4)+fabs(a[2])*Ltot))/Ptot);
    y = ((a[0] + (a[1]-a[0])*temp)*Ptot)/a.back() ;
    dyda[0]= (1 - temp)*Ptot/a.back();
    dyda[1]= temp*Ptot/a.back();
    dyda[2]=Ptot*(a[1]-a[0])*(0.5-((((a[2]+Ptot-Ltot)/2)+Ltot) /
                              sqrt((((a[2]+Ptot-Ltot)*(a[2]+Ptot-Ltot))/4)+Ltot*a[2])))/(Ptot*a.back());
}

double fit_one_site_dilution(const double Ltot, const double Ptot, const std::vector<double> &a)
{
    double temp = ((Ltot/Ptot)+(((Ptot-Ltot+fabs(a[2]))/2) -
                   sqrt((((Ptot-Ltot+fabs(a[2]))*(Ptot-Ltot+fabs(a[2])))/4)+fabs(a[2])*Ltot))/Ptot);
    return ((a[0] + (a[1]-a[0])*temp)*Ptot)/a.back() ;
}

//Two-site
// -------
// a0-a2: P - PL2 signals
// a3-a4: Kd1 - Kd2
// a5: [L]start
// a6: [P]start

double ligand_free_two_sites_dilution(const double start_ligand, const double protein_conc, const std::vector<double> &a, const double Lfree);
double protein_free_two_sites_dilution(const double L, const double x2, const std::vector<double> &a);
double secord_deriv(double ligand, unsigned int idx, double protein, VecDoub_I &a,
                                      double fitfunc(double, double, VecDoub_I &));

void fit_two_sites_dilution(const double Ltot, const double Ptot, const std::vector<double> &a, double &y, std::vector<double> &dyda)
{
    double root{}, start{a[a.size()-2]}; //a[5] is start conc of the datasets/conc_vector, since we need it to make zbrent

    try {
        root=zbrent(ligand_free_two_sites_dilution,start,Ltot,1e-12, Ltot,Ptot,a);
    } catch (int) {
    }
    double P=protein_free_two_sites_dilution(root,Ptot,a);
    double PL=P*root/fabs(a[3]);
    double PL2=PL*root/fabs(a[4]);
    y = ((a[0]*P + a[1]*PL + a[2]*PL2)/Ptot)*Ptot/a.back();
    dyda[0] = P/a.back();
    dyda[1] = PL/a.back();
    dyda[2] = PL2/a.back();
    for (unsigned int i=3;i<dyda.size();i++)
        dyda[i]=secord_deriv(Ltot, i, Ptot, a, fit_two_sites_dilution);
}

double fit_two_sites_dilution(const double Ltot, const double Ptot, const std::vector<double> &a)
{
    double root{}, Lstart{a[a.size()-2]}; //a[5] is start conc of the datasets/conc_vector, since we need it to make zbrent
    try {
        root=zbrent(ligand_free_two_sites_dilution, Lstart, Ltot, 1e-12, Ltot, Ptot,a);
    } catch (int) {

    }
    double P = protein_free_two_sites_dilution(root,Ptot,a);
    double PL = P*root/fabs(a[3]);
    double PL2 = PL*root/fabs(a[4]);
    return (a[0]*P + a[1]*PL + a[2] * PL2)/a.back();
}

double ligand_free_two_sites_dilution(const double Ltot, const double Ptot, const std::vector<double> &a, const double L) {
    double sum = L/fabs(a[3]) + 2*L*L/fabs(a[3]*a[4]);
    return Ltot - L - Ptot*sum/(1.+sum);
}

double protein_free_two_sites_dilution(const double L, const double Ptot, const std::vector<double> &a){
    return Ptot/(1+(L/fabs(a[3]))+(L*L)/fabs(a[3]*a[4]));
}


// Four-site
// ---------
// a0-a4: P - PL4 signals
// a5-a8: Kd1 - Kd4
// a9: [L]start
// a10: [P]start


double ligand_free_four_sites_dilution(double Ltot, double Ptot, VecDoub_I &a, const double L);
double protein_free_four_sites_dilution(double L, double Ptot, VecDoub_I &a);
double secord_deriv(double ligand, unsigned int idx, double protein, VecDoub_I &a,
                                      double fitfunc(double, double, VecDoub_I &));

double fit_four_sites_dilution(const double Ltot, const double Ptot, const std::vector<double> &a)
{
    double root{}, Lstart{a[a.size()-2]};
    try {
        root=zbrent(ligand_free_four_sites_dilution, Lstart, Ltot, 1e-12, Ltot, Ptot, a);
    } catch (int) {

    }
    double P = protein_free_four_sites_dilution(root, Ptot, a);
    double PL=P*root/fabs(a[5]);
    double PL2=PL*root/fabs(a[6]);
    double PL3=PL2*root/fabs(a[7]);
    double PL4=PL3*root/fabs(a[8]);
    return ((a[0]*P + a[1]*PL + a[2]*PL2 + a[3]*PL3 + a[4]*PL4)/Ptot)*Ptot/a.back();
}

void fit_four_sites_dilution(const double Ltot, const double Ptot, const std::vector<double> &a, double &y, std::vector<double> &dyda)
{
    double root{}, Lstart{a[a.size()-2]};

    try {
        root=zbrent(ligand_free_four_sites_dilution, Lstart, Ltot, 1e-12, Ltot, Ptot, a);
    } catch (int) {
    }
    double P=protein_free_four_sites_dilution(root,Ptot,a);
    double PL=P*root/fabs(a[5]);
    double PL2=PL*root/fabs(a[6]);
    double PL3=PL2*root/fabs(a[7]);
    double PL4=PL3*root/fabs(a[8]);
    y = (a[0]*P + a[1]*PL + a[2]*PL2 + a[3]*PL3 + a[4]*PL4)/a.back();

    dyda[0] = P/a.back();
    dyda[1] = PL/a.back();
    dyda[2] = PL2/a.back();
    dyda[3] = PL3/a.back();
    dyda[4] = PL4/a.back();
    for (unsigned int i=5;i<dyda.size();i++)
        dyda[i]=secord_deriv(Ltot, i, Ptot, a, fit_four_sites_dilution);
}

double ligand_free_four_sites_dilution(double Ltot, double Ptot, VecDoub_I &a, double L){
    double sum = L/fabs(a[5]) + 2*L*L/fabs(a[5]*a[6]) + 3*L*L*L/fabs(a[5]*a[6]*a[7]) +
           4*L*L*L*L/fabs(a[5]*a[6]*a[7]*a[8]);
    return Ltot - L - Ptot*sum/(1.+sum);
}

double protein_free_four_sites_dilution(double L, double Ptot, VecDoub_I &a) {
    return Ptot/(1 + L/fabs(a[5]) + L*L/fabs(a[5]*a[6]) + L*L*L/fabs(a[5]*a[6]*a[7]) +
            L*L*L*L/fabs(a[5]*a[6]*a[7]*a[8]));
}



// Two site binding + Competitive Binding
// --------------------------------------
// a0-a1: Comp - CompL signals
// a2-a4: Kd1, Kd2, KdComp
// a5: [L]start
// a6: [P]start
// a7: [Comp]start

double secord_deriv_comp(const double Ltot, unsigned int idx, const double Ptot, const std::vector<double> &a,
                           const double comp);

void fit_comp_dilution(const double Ltot, const double Ptot, const double CompTot, const std::vector<double> &a, double &y, std::vector<double> &dyda)
{
    double I{},IL{}, root{}, Lstart{a[a.size()-3]}; //a[5] is start conc of the datasets/conc_vector, since we need it to make zbrent
    try {
        root=zbrent2(ligand_free_comp,Lstart,Ltot,1e-12, Ltot,Ptot,CompTot,a);
    } catch (int) {

    }
    I = (fabs(a[4])*CompTot)/(fabs(a[4])+root);
    IL = CompTot-I;
    y = (a[0]*I + a[1]*IL)/a.back();

    for (unsigned int i=0;i<dyda.size();i++)
        dyda[i]=secord_deriv_comp(Ltot, i, Ptot, a, CompTot);
}

double fit_comp_dilution(const double Ltot, const double Ptot, const double CompTot, const std::vector<double> &a)
{
    double I{},IL{}, root{}, Lstart{a[a.size()-3]};
    try {
        root=zbrent2(ligand_free_comp, Lstart, Ltot, 1e-12, Ltot,Ptot, CompTot, a);
    } catch (int) {

    }
    I=(fabs(a[4])*CompTot)/(fabs(a[4])+root);
    IL=CompTot - I;
    return (a[0]*I + a[1]*IL)/a.back();
}

double ligand_free_comp(const double start_ligand, const double protein_conc, const double comp_conc, const std::vector<double> &a, const double Lfree)
{
    return Lfree+comp_conc*Lfree/(fabs(a[4])+Lfree) + (protein_conc/(1+Lfree/fabs(a[2])+2*(Lfree*Lfree)/fabs(a[2]*a[3]))) *
            (Lfree/fabs(a[2])+2*(Lfree*Lfree)/fabs(a[2]*a[3])) - start_ligand;
}

//*********************************************************************************************
// THE FIRST FUNCTIONS BELOW ARE USED TO FIT DATA WHERE THE SIGNAL DOES NOT DEPEND ON CONCENTRATION
// EXAMPLES: CHEMICAL SHIFTS, FLUORESCENCE ANISOTROPY
//*********************************************************************************************

//x1=ligand concentration, x2=proteinconcetration, y=responce, a vector yp(fritt protein) är a[o], y_pl=a[1](signal från PL) och kd= a[2]
// One-site (dilution not important for signal)
// --------------------------------------------
// a0-a1: P - PL signals
// a2: Kd1 - Kd2
// a3: [L]start
// a4: [P]start

void fit_one_site(const double Ltot, const double Ptot, const std::vector<double> &a, double &y, std::vector<double> &dyda)
{
    double temp = ((Ltot/Ptot)+(((Ptot-Ltot+fabs(a[2]))/2) -
          sqrt((((Ptot-Ltot+fabs(a[2]))*(Ptot-Ltot+fabs(a[2])))/4)+fabs(a[2])*Ltot))/Ptot);

    y = a[0] + (a[1]-a[0])*temp;

    dyda[0] = 1 - temp;
    dyda[1] = temp;
    dyda[2] = (a[1]-a[0])*(0.5-((((a[2]+Ptot-Ltot)/2)+Ptot) /
                         sqrt((((a[2]+Ptot-Ltot)*(a[2]+Ptot-Ltot))/4)+Ltot*a[2])))/Ptot;
}

double fit_one_site(const double Ltot, const double Ptot, const std::vector<double> &a)
{
    double temp=((Ltot/Ptot)+(((Ptot-Ltot+a[2])/2) -
          sqrt((((Ptot-Ltot+a[2])*(Ptot-Ltot+a[2]))/4)+a[2]*Ltot)) /Ptot);
    return a[0] + (a[1]-a[0])*temp;
}

// Two-site (dilution not important for signal)
// --------------------------------------------
// a0-a2: P - PL2 signals
// a3-a4: Kd1 - Kd2
// a5: [L]start
// a6: [P]start

double ligand_free_two_sites(const double start_ligand, const double protein_conc, const std::vector<double> &a, const double Lfree);
double protein_free_two_sites(const double L, const double x2, const std::vector<double> &a);
double secord_deriv(double ligand, unsigned int idx, double protein, VecDoub_I &a,
                                      double fitfunc(double, double, VecDoub_I &));
double fit_two_sites(const double Ltot, const double Ptot, const std::vector<double> &a)
{
    double P{}, root{}, y{}, Lstart{a[a.size()-2]}; //a[5] is start conc of the datasets/conc_vector, since we need it to make zbrent
    try {
        root=zbrent(ligand_free_two_sites, Lstart, Ltot, 1e-12, Ltot, Ptot, a);
    } catch (int) {

    }
    P=protein_free_two_sites(root, Ptot, a);
    double PL{}, PL2{};
    PL=P*root/fabs(a[3]);
    PL2=PL*root/fabs(a[4]);
    y = (a[0]*P + a[1]*PL + a[2] * PL2)/Ptot;
    return y;
}

void fit_two_sites(const double ligand_conc, const double Ptot, const std::vector<double> &a, double &y, std::vector<double> &dyda)
{
    double P{}, root{}, start{a[a.size()-2]}; //a[5] is start conc of the datasets/conc_vector, since we need it to make zbrent

    try {
        root=zbrent(ligand_free_two_sites, start, ligand_conc, 1e-12, ligand_conc, Ptot, a);
    } catch (int) {

    }
    P=protein_free_two_sites(root, Ptot, a);
    double PL{}, PL2{};
    PL=P*root/fabs(a[3]);
    PL2=PL*root/fabs(a[4]);
    y = (a[0]*P + a[1]*PL + a[2] * PL2)/Ptot;
    dyda[0] = P/Ptot;
    dyda[1] = PL/Ptot;
    dyda[2] = PL2/Ptot;
    for (unsigned int i=3;i<dyda.size();i++)
        dyda[i]=secord_deriv(ligand_conc, i, Ptot, a, fit_two_sites);
}

double ligand_free_two_sites(const double start_ligand, const double protein_conc, const std::vector<double> &a, const double Lfree)
{
    return (Lfree*Lfree*Lfree/fabs(a[3]*a[4]))+(Lfree*Lfree)*((1/fabs(a[3]))+(2*protein_conc/(a[3]*a[4]))-1/fabs(a[3]*a[4]))+
            Lfree*(1+(protein_conc/fabs(a[3]))-(start_ligand/fabs(a[3])))-start_ligand;
}

double protein_free_two_sites(const double L, const double x2, const std::vector<double> &a)
{
    return x2/(1+(L/a[3])+(L*L)/(a[3]*a[4]));
}

double fit_four_sites(const double Ltot, const double Ptot, const std::vector<double> &a)
{
    double root{}, Lstart{a[a.size()-2]};
    try {
        root=zbrent(ligand_free_four_sites_dilution, Lstart, Ltot, 1e-12, Ltot, Ptot, a);
    } catch (int) {

    }
    double P = protein_free_four_sites_dilution(root, Ptot, a);
    double PL=P*root/fabs(a[5]);
    double PL2=PL*root/fabs(a[6]);
    double PL3=PL2*root/fabs(a[7]);
    double PL4=PL3*root/fabs(a[8]);
    return (a[0]*P + a[1]*PL + a[2]*PL2 + a[3]*PL3 + a[4]*PL4)/Ptot;
}

void fit_four_sites(const double Ltot, const double Ptot, const std::vector<double> &a, double &y, std::vector<double> &dyda)
{
    double root{}, start{a[a.size()-2]};

    try {
        root=zbrent(ligand_free_four_sites_dilution,start,Ltot,1e-12, Ltot,Ptot,a);
    } catch (int) {
    }
    double P=protein_free_four_sites_dilution(root,Ptot,a);
    double PL=P*root/fabs(a[5]);
    double PL2=PL*root/fabs(a[6]);
    double PL3=PL2*root/fabs(a[7]);
    double PL4=PL3*root/fabs(a[8]);
    y = (a[0]*P + a[1]*PL + a[2]*PL2 + a[3]*PL3 + a[4]*PL4)/Ptot;

    dyda[0] = P/Ptot;
    dyda[1] = PL/Ptot;
    dyda[2] = PL2/Ptot;
    dyda[3] = PL3/Ptot;
    dyda[4] = PL4/Ptot;
    for (unsigned int i=5;i<dyda.size();i++)
        dyda[i]=secord_deriv(Ltot, i, Ptot, a, fit_four_sites);
}


void fit_comp(const double Ltot, const double Ptot, const double CompTot, const std::vector<double> &a, double &y, std::vector<double> &dyda)
{
    double I{},IL{}, root{}, Lstart{a[a.size()-3]}; //a[5] is start conc of the datasets/conc_vector, since we need it to make zbrent
    try {
        root=zbrent2(ligand_free_comp,Lstart,Ltot,1e-12, Ltot,Ptot,CompTot,a);
    } catch (int) {

    }
    I = (fabs(a[4])*CompTot)/(fabs(a[4])+root);
    IL = CompTot-I;
    y = (a[0]*I + a[1]*IL)/CompTot;

    for (unsigned int i=0;i<dyda.size();i++)
        dyda[i]=secord_deriv_comp(Ltot, i, Ptot, a, CompTot);
}

double fit_comp(const double Ltot, const double Ptot, const double CompTot, const std::vector<double> &a)
{
    double I{},IL{}, root{}, Lstart{a[a.size()-3]};
    try {
        root=zbrent2(ligand_free_comp, Lstart, Ltot, 1e-12, Ltot,Ptot, CompTot, a);
    } catch (int) {

    }
    I=(fabs(a[4])*CompTot)/(fabs(a[4])+root);
    IL=CompTot - I;
    return (a[0]*I + a[1]*IL)/CompTot;
}

double secord_deriv_comp(const double ligand, unsigned int idx, const double protein, const std::vector<double> &a, const double comp)
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
    yh = fit_comp_dilution(ligand,protein,comp, a_temp);
    a_temp[idx] = aidx - h;
    ymh = fit_comp_dilution(ligand,protein,comp, a_temp);
    double diff = yh - ymh;
    a_temp[idx] = aidx;
    return diff/(2*h);
}

double secord_deriv(double Ltot, unsigned int idx, double Ptot, VecDoub_I &a,
                                      double fitfunc(double, double, VecDoub_I &))
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
    yh = fitfunc(Ltot, Ptot, a_temp);
    a_temp[idx] = aidx - h;
    ymh = fitfunc(Ltot, Ptot, a_temp);
    double diff = yh - ymh;
    a_temp[idx] = aidx;
    return diff/(2*h);
}

bool gaussj(MatDoub_IO &a, MatDoub_IO &b)
{
    int i,icol,irow,j,k,l,ll,n=a.nrows(),m=b.ncols();
    double big,dum,pivinv;
    VecInt indxc(n),indxr(n),ipiv(n);
    for (j=0;j<n;j++) ipiv[j]=0;
    for (i=0;i<n;i++) {
        big=0.0;
        for (j=0;j<n;j++)
            if (ipiv[j] != 1)
                for (k=0;k<n;k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big=fabs(a[j][k]);
                            irow=j;
                            icol=k;
                        }
                    }
                }
        ++(ipiv[icol]);
        if (irow != icol) {
            for (l=0;l<n;l++) std::swap(a[irow][l],a[icol][l]);
            for (l=0;l<m;l++) std::swap(b[irow][l],b[icol][l]);
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
                std::swap(a[k][indxr[l]],a[k][indxc[l]]);
    }
    return true;
}


template<class T>
inline T SIGN(const T &a, const T &b)
    {return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

template <class T>
double zbrent(T &func, const double x1, const double x2, const double tol, const double xx1, const double xx2, const std::vector<double> &aa)
{
    const int ITMAX=100;
    const double EPS=std::numeric_limits<double>::epsilon();
    double a=x1,b=x2,c=x2,d,e,fa=func(xx1,xx2, aa, a),fb=func(xx1,xx2, aa, b),fc,p,q,r,s,tol1,xm;
    //std::cout << "fa = " << fa << ", fb = " << fb <<  ", x1 = "<<x1 << ", x2 = " << x2 << std::endl;
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        throw("Root must be bracketed in zbrent");
    fc=fb;
    for (int iter=0;iter<ITMAX;iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*EPS*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) return b;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
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
            p=fabs(p);
            double min1=3.0*xm*q-fabs(tol1*q);
            double min2=fabs(e*q);
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
        if (fabs(d) > tol1)
            b += d;
        else
            b += SIGN(tol1,xm);
            fb=func(xx1,xx2, aa, b);
    }
    throw("Maximum number of iterations exceeded in zbrent");
}

template <class T>
double zbrent2(T &func, const double x1, const double x2, const double tol, const double xx1, const double xx2, const double xx3,const std::vector<double> &aa)
{
    const int ITMAX=100;
    const double EPS=std::numeric_limits<double>::epsilon();
    double a=x1,b=x2,c=x2,d,e,fa=func(xx1,xx2,xx3, aa, a),fb=func(xx1,xx2,xx3, aa, b),fc,p,q,r,s,tol1,xm;
    //std::cout << "fa = " << fa << ", fb = " << fb <<  ", x1 = "<<x1 << ", x2 = " << x2 << std::endl;
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        throw("Root must be bracketed in zbrent");
    fc=fb;
    for (int iter=0;iter<ITMAX;iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*EPS*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) return b;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
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
            p=fabs(p);
            double min1=3.0*xm*q-fabs(tol1*q);
            double min2=fabs(e*q);
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
        if (fabs(d) > tol1)
            b += d;
        else
            b += SIGN(tol1,xm);
            fb=func(xx1,xx2, xx3,aa, b);
    }
    throw("Maximum number of iterations exceeded in zbrent");
}

//All guess parameter functions is when the software estimates a starting guess
//All guess parameters_range is when the user has guessed within a range.

void guess_parameters_one_site(Data *&data, QVector<double> kd_testVector, QString mode)
{
    double y_calc{};
    double chi2{-1.0};
    double chi2_temp{};
    int best_index_kd{}, best_index_responce{};
    QVector<double> test_last_responce{};
    double range{data->responceVector.last()/10};
    for (int i=0;i<11;i++)
        test_last_responce.push_back(data->responceVector.last()+i*range);

    double (*dblptr)(const double, const double, const std::vector<double> &){fit_one_site_dilution};
    if (mode=="chemshift")
        dblptr = fit_one_site;
    for (int i=0; i<kd_testVector.size();i++){
        std::vector<double> a_temp{data->a};
        a_temp[2]=kd_testVector[i];
        for (int h=0;h<test_last_responce.size();h++){
            a_temp[1]=test_last_responce[h];
            for (int i=0; i<data->responceVector.size();i++){
                y_calc=dblptr(data->concVector[i],data->protein_conc_vector[i], a_temp);
                chi2_temp += ((y_calc - data->responceVector[i]) *
                              (y_calc - data->responceVector[i]))
                        / (data->error_vector[i]*data->error_vector[i]);
            }

            if(!std::isnan(chi2_temp)){
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

void guess_parameters_range(Data *&data, QVector<double> &kd, QVector<double> &start, QVector<double> &end, QString mode)
{
    double y_calc{};
    double chi2{};
    double chi2_temp{};
    int best_index_kd{}, best_index_start_responce{}, best_index_end_responce{};
    std::vector<double> a_temp{double(),double(),double()};
    double (*dblptr)(const double, const double, const std::vector<double> &){fit_one_site_dilution};
    if (mode=="chemshift")
        dblptr = fit_one_site;
    for (int i = 0; i < kd.size();i++) {
        a_temp[2]=kd[i];
        for (int j=0; j < start.size();j++){
            a_temp[0]=start[j];
            for (int l=0;l<end.size();l++){
                a_temp[1]=end[l];
                for (int k=0; k<data->responceVector.size();k++){
                    y_calc=dblptr(data->concVector[k],data->protein_conc_vector[k], a_temp);
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

void guess_parameters_range_2_sites(Data *&data, QVector<double> &kd, QVector<double> &kd2, QVector<double> &start, QVector<double> &inter, QVector<double> &end, QString mode)
{
    double y_calc{};
    double chi2{-1.0};
    double chi2_temp{};
    int best_index_kd{}, best_index_start_responce{}, best_index_end_responce{}, best_index_inter_responce{}, best_index_kd2{};
    std::vector<double> a_temp{data->a};
    double (*dblptr)(const double, const double, const std::vector<double> &){fit_two_sites_dilution};
    if (mode=="chemshift")
        dblptr = fit_two_sites;
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
                            y_calc=dblptr(data->concVector[k],data->protein_conc_vector[k], a_temp);
                            chi2_temp += ((y_calc - data->responceVector[k]) *
                                          (y_calc - data->responceVector[k]))
                                    / (data->error_vector[k]*data->error_vector[k]);
                        }
                        if(!std::isnan(chi2_temp)){
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

void guess_parameters_two_sites(Data *&data, QVector<double> kd_testVector, QProgressDialog &progress, QString mode)
{
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
    double (*dblptr)(const double, const double, const std::vector<double> &){fit_two_sites_dilution};
    if (mode=="chemshift")
        dblptr = fit_two_sites;
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
                        y_calc=dblptr(data->concVector[k],data->protein_conc_vector[k], a_temp);
                        chi2_temp += ((y_calc - data->responceVector[k]) *
                                      (y_calc - data->responceVector[k]))
                                / (data->error_vector[k]*data->error_vector[k]);
                    }
                    if(!std::isnan(chi2_temp)){
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

void guess_parameters_range_four_sites(Data *&data, QVector<double> &kd1_testVector, QVector<double> &kd2_testVector, QVector<double> &kd3_testVector,
                                       QVector<double> &kd4_testVector, QVector<double> &p_resp_testVector, QVector<double> &pl_resp_testVector,
                                       QVector<double> &pl2_resp_testVector, QVector<double> &pl3_resp_testVector, QVector<double> &pl4_resp_testVector,
                                       QString mode)
{
    double chi2 = std::numeric_limits<double>::max();
    VecInt bestidx{};
    double (*dblptr)(const double, const double, const std::vector<double> &){fit_four_sites_dilution};
    if (mode=="chemshift")
        dblptr = fit_four_sites;
    for (int i0=0; i0<p_resp_testVector.size(); ++i0)
        for (int i1=0; i1<pl_resp_testVector.size(); ++i1)
            for (int i2=0; i2<pl2_resp_testVector.size(); ++i2)
                for (int i3=0; i3<pl3_resp_testVector.size(); ++i3)
                    for (int i4=0; i4<pl4_resp_testVector.size(); ++i4)
                        for (int i5=0; i5<kd1_testVector.size(); ++i5)
                            for (int i6=0; i6<kd2_testVector.size(); ++i6)
                                for (int i7=0; i7<kd3_testVector.size(); ++i7)
                                    for (int i8=0; i8<kd4_testVector.size(); ++i8) {
                                        double chi2_temp = 0.;
                                        VecDoub a_temp = {p_resp_testVector[i0], pl_resp_testVector[i1],
                                                  pl2_resp_testVector[i2], pl3_resp_testVector[i3],
                                                  pl4_resp_testVector[i4], kd1_testVector[i5],
                                                  kd2_testVector[i6], kd3_testVector[i7],
                                                  kd4_testVector[i8], data->a[9], data->a[10]};
                                        for (int j=0; j<data->concVector.size(); ++j) {
                                             double y = dblptr(data->concVector[j],
                                                               data->protein_conc_vector[j], a_temp);
                                            chi2_temp += ((y - data->responceVector[j]) *
                                                          (y - data->responceVector[j]))
                                                    / (data->error_vector[j]*data->error_vector[j]);
                                        }
                                        if (chi2_temp<chi2) {
                                            bestidx = {i0, i1, i2, i3, i4, i5, i6, i7, i8};
                                            chi2 = chi2_temp;
                                        }
                                    }
    data->a[0]=p_resp_testVector[bestidx[0]];
    data->a[1]=pl_resp_testVector[bestidx[1]];
    data->a[2]=pl2_resp_testVector[bestidx[2]];
    data->a[3]=pl3_resp_testVector[bestidx[3]];
    data->a[4]=pl4_resp_testVector[bestidx[4]];
    data->a[5]=kd1_testVector[bestidx[5]];
    data->a[6]=kd2_testVector[bestidx[6]];
    data->a[7]=kd3_testVector[bestidx[7]];
    data->a[8]=kd4_testVector[bestidx[8]];
}


// Fix progress later
void guess_parameters_four_sites(Data *&data, QProgressDialog &progress, QString mode)
{
    double chi2 = std::numeric_limits<double>::max();
    VecInt bestidx{};
    double (*dblptr)(const double, const double, const std::vector<double> &){fit_four_sites_dilution};
    if (mode=="chemshift")
        dblptr = fit_four_sites;

    int maxi1{10}, maxi2{maxi1}, maxi3{maxi1}, maxi4{maxi1}, maxi5{6}, maxi6{maxi5}, maxi7{maxi5}, maxi8{maxi5};
    VecDoub response_range{10};
    for (int i=0; i<2; ++i) response_range[i] = data->responceVector[0]+i*(data->responceVector.back()-data->responceVector[0])/(maxi1-1);
    VecDoub kd_range{0.001, 0.01, 0.1, 1., 10., 100.};
    for (int i1=0; i1<maxi1; ++i1)for (int i1=0; i1<maxi1; ++i1) {
        progress.setValue(i1+1);
        for (int i2=0; i2<maxi2; ++i2)
            for (int i3=0; i3<maxi3; ++i3)
                for (int i4=0; i4<maxi4; ++i4)
                    for (int i5=0; i5<maxi5; ++i5)
                        for (int i6=0; i6<maxi6; ++i6)
                            for (int i7=0; i7<maxi7; ++i7)
                                for (int i8=0; i8<maxi8; ++i8) {
                                    double chi2_temp = 0.;
                                    VecDoub a_temp = {data->concVector[0], response_range[i1],
                                                      response_range[i2], response_range[i3],
                                                      response_range[i4], kd_range[i5],
                                                      kd_range[i6], kd_range[i7],
                                                      kd_range[i8], data->a[9], data->a[10]};
                                    for (int j=0; j<data->concVector.size(); ++j) {
                                        double y = dblptr(data->concVector[j],
                                                          data->protein_conc_vector[j], a_temp);
                                        chi2_temp += ((y - data->responceVector[j]) *
                                                      (y - data->responceVector[j]))
                                                / (data->error_vector[j]*data->error_vector[j]);
                                    }
                                    if (chi2_temp<chi2) {
                                        bestidx = {i1, i2, i3, i4, i5, i6, i7, i8};
                                        chi2 = chi2_temp;
                                    }
                                }
    }
    progress.setValue(maxi1+1);
    data->a[0]=data->concVector[0];
    data->a[1]=response_range[bestidx[0]];
    data->a[2]=response_range[bestidx[1]];
    data->a[3]=response_range[bestidx[2]];
    data->a[4]=response_range[bestidx[3]];
    data->a[5]=kd_range[bestidx[4]];
    data->a[6]=kd_range[bestidx[5]];
    data->a[7]=kd_range[bestidx[6]];
    data->a[8]=kd_range[bestidx[7]];
}


void guess_parameters_comp(Data *&data, QVector<double> kd_testVector, QProgressDialog &progress, QString mode)
{
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
    double (*dblptr)(const double, const double, const double, const std::vector<double> &){fit_comp_dilution};
    if (mode=="chemshift")
        dblptr = fit_comp;
    for (int g = 0; g < test_last_responce.size();g++) {
        a_temp[1]=test_last_responce[g];
        progress.setValue(g+1);
        for (int i=0; i<kd_testVector.size();i++){
            a_temp[2]=kd_testVector[i];
            for (int j=i;j<kd_testVector.size();j++){
                a_temp[3]=kd_testVector[j];
                for (int c=0;c<kd_testVector.size();c++){
                    a_temp[4]=kd_testVector[c];
                    for (int k=0; k<data->responceVector.size();k++) {
                        y_calc=dblptr(data->concVector[k],data->protein_conc_vector[k],data->comp_vector[k], a_temp);
                        chi2_temp += ((y_calc - data->responceVector[k]) *
                                      (y_calc - data->responceVector[k]))
                                / (data->error_vector[k]*data->error_vector[k]);
                    }
                    if(!std::isnan(chi2_temp)){
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

void guess_parameters_comp_range(Data *&data, QVector<double> &kd, QVector<double> &kd2, QVector<double> &kdc, QVector<double> &start, QVector<double> &end, QString mode)
{
    double y_calc{};
    double chi2{-1.0};
    double chi2_temp{};
    int best_index_kd1{}, best_index_start_responce{}, best_index_end_responce{}, best_index_kdc{}, best_index_kd2{};
    std::vector<double> a_temp{data->a};
    double (*dblptr)(const double, const double, const double, const std::vector<double> &){fit_comp_dilution};
    if (mode=="chemshift")
        dblptr = fit_comp;
    for (int i = 0; i < kd.size();i++) {
        a_temp[2]=kd[i];
        for (int j=0; j < start.size();j++) {
            a_temp[0]=start[j];
            for (int l=0; l < end.size();l++) {
                a_temp[1]=end[l];
                for (int m=0;m<kdc.size();m++) {
                    a_temp[4]=kdc[m];
                    for (int n=0;n<kd2.size();n++) {
                        a_temp[3]=kd2[n];
                        for (int k=0; k<data->responceVector.size();k++){
                            y_calc=dblptr(data->concVector[k],data->protein_conc_vector[k],data->comp_vector[k], a_temp);
                            chi2_temp += ((y_calc - data->responceVector[k]) *
                                          (y_calc - data->responceVector[k]))
                                    / (data->error_vector[k]*data->error_vector[k]);
                        }
                        if(!std::isnan(chi2_temp)) {
                            if (chi2<0) {
                                chi2=chi2_temp;
                                best_index_end_responce=l;
                                best_index_start_responce=j;
                                best_index_kd1=i;
                                best_index_kdc=m;
                                best_index_kd2=n;
                            }
                            else if(chi2>chi2_temp) {
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

void guess_parameters_cpmg(Data *&data, QProgressDialog &progress)
{
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
                    if(!std::isnan(chi2_temp)){
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


void guess_parameters_cpmg_range(Data *&data, QVector<double> &r20, QVector<double> &kex, QVector<double> &pb, QVector<double> &dw)
{
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
                    if(!std::isnan(chi2_temp)){
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

/*
double secord_deriv(const double ligand, unsigned int idx, const double protein, const std::vector<double> &a)
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
    double yh = fit_two_sites(ligand,protein, a_temp);
    a_temp[idx] = aidx - h;
    double ymh = fit_two_sites(ligand,protein, a_temp);
    double diff = yh - ymh;
    a_temp[idx] = aidx;
    return diff/(2*h);
}
*/

template<class T>
inline T SQR(const T a) {return a*a;}

double RexCR(double R20, double DW, double PB, double KEX, double nuCP)
{
// R20 = (R2A + R2B)/2
// Note that order of arguments does not correspond to order of a-vector!!
// Also note that 2*tauCP is separation of 180-pulses

    std::complex<double> psi, zeta, np, nm, Dp, Dm;
    double tmp;
    double PA = 1.-PB;
    double tauCP = 0.25/nuCP;

    psi = SQR(-PA*KEX+PB*KEX)-SQR(DW)+4.*PA*PB*SQR(KEX);
    zeta = 2.*DW*(-PA*KEX+PB*KEX);
    np = sqrt(2.)*tauCP*sqrt(sqrt(SQR(psi) + SQR(zeta)) + psi);
    nm = sqrt(2.)*tauCP*sqrt(sqrt(SQR(psi) + SQR(zeta)) - psi);
    Dp = 0.5*((psi + 2.*SQR(DW))/sqrt(SQR(psi) + SQR(zeta)) + 1.);
    Dm = 0.5*((psi + 2.*SQR(DW))/sqrt(SQR(psi) + SQR(zeta)) - 1.);
    tmp = fabs(real(acosh_complex(Dp*cosh(np)-Dm*cos(nm))));
    return 0.5*(2.*R20+KEX-1./(2.*tauCP)*tmp);
}

void carverrichards(double x, VecDoub_I &a, double &y, VecDoub_O &dyda)
{
    const double STEP = 1.0e-6;
    double R20 = a[0];
    double KEX = a[1];
    double PB = a[2];
    double DW = a[3];

    y = RexCR(R20, DW, PB, KEX, x);

    dyda[0] = 1.0;
    dyda[1] = (RexCR(R20, DW, PB, (1.+STEP)*KEX, x) - y ) / (STEP*KEX);
    dyda[2] = (RexCR(R20, DW, (1.+STEP)*PB, KEX, x) - y ) / (STEP*PB);
    dyda[3] = (RexCR(R20, (1.+STEP)*DW, PB, KEX, x) - y ) / (STEP*DW);
}

double carverrichards(double x, VecDoub_I &a)
{
    double R20 = a[0];
    double KEX = a[1];
    double PB = a[2];
    double DW = a[3];
    return RexCR(R20, DW, PB, KEX, x);
}

inline std::complex<double> acosh_complex(std::complex<double> z) {
    return log(z + sqrt(z-1.)*sqrt(z+1.));
}

double CPMG_kd(VecDoub_I &a){
    double PB = a[2];
    double Kb = a[1]-PB*a[1];
    double Kf_L = a[1]-Kb;
    return (Kb * a.back())/ Kf_L;
}
