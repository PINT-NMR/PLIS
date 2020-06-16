#ifndef FITMRQ_H
#define FITMRQ_H

#include "calculations.h"
#include "nr3.h"

// 200530
// - added constructors that take an ia vector

struct Fitmrq {
    static const int NDONE=4, ITMAX=1000;
    int ndat, ma, mfit;
    VecDoub_I x, x2, y, volume, sig;
    double tol;
    void (*funcs)(const double, const double, VecDoub_I &, double &, VecDoub_O &);
    double (*funcs2)(const double, const double, const std::vector<double> &);
	VecBool ia;
    MatDoub alpha;
    VecDoub a;
    MatDoub covar;
    double chisq;
    bool did_gauss_work;

    Fitmrq(VecDoub_I &xx, VecDoub_I &xx2,VecDoub_I &yy, VecDoub_I & v,VecDoub_I &ssig, VecDoub_I &aa,
    void funks(const double, const double, VecDoub_I &, double &, VecDoub_O &), double funks2(const double, const double, const std::vector<double> &),
           const double TOL=1.e-3) : ndat(xx.size()), ma(aa.size()), x(xx), x2(xx2),y(yy), volume(v), sig(ssig),
    tol(TOL), funcs(funks), funcs2(funks2),ia(ma), alpha(ma,ma), a(aa), covar(ma,ma) {
        for (int i=0;i<ma;i++)
            ia[i] = true;
    }

    Fitmrq(VecDoub_I &xx, VecDoub_I &xx2,VecDoub_I &yy, VecDoub_I & v,VecDoub_I &ssig, VecDoub_I &aa, VecBool &iia,
    void funks(const double, const double, VecDoub_I &, double &, VecDoub_O &), double funks2(const double, const double, const std::vector<double> &),
           const double TOL=1.e-3) : ndat(xx.size()), ma(aa.size()), x(xx), x2(xx2),y(yy), volume(v), sig(ssig),
    tol(TOL), funcs(funks), funcs2(funks2),ia(ma), alpha(ma,ma), a(aa), covar(ma,ma) {
        if (iia.size() == ia.size())
            ia = iia;
        else
            for (int i=0;i<ma;i++)
                ia[i] = true;
    }

    int fit() {
        int j,k,l,iter,done=0;
        double alamda=.001,ochisq;
		VecDoub atry(ma),beta(ma),da(ma);
		mfit=0;
		for (j=0;j<ma;j++) if (ia[j]) mfit++;
		MatDoub oneda(mfit,1), temp(mfit,mfit);
		mrqcof(a,alpha,beta);
		for (j=0;j<ma;j++) atry[j]=a[j];
		ochisq=chisq;
		for (iter=0;iter<ITMAX;iter++) {
			if (done==NDONE) alamda=0.;
			for (j=0;j<mfit;j++) {
				for (k=0;k<mfit;k++) covar[j][k]=alpha[j][k];
				covar[j][j]=alpha[j][j]*(1.0+alamda);
				for (k=0;k<mfit;k++) temp[j][k]=covar[j][k];
				oneda[j][0]=beta[j];
			}
            did_gauss_work=gaussj(temp,oneda);
            if (did_gauss_work==false)
                return 1;
			for (j=0;j<mfit;j++) {
				for (k=0;k<mfit;k++) covar[j][k]=temp[j][k];
				da[j]=oneda[j][0];
			}
			if (done==NDONE) {
				covsrt(covar);
				covsrt(alpha);
                return 0;
			}
			for (j=0,l=0;l<ma;l++)
				if (ia[l]) atry[l]=a[l]+da[j++];
			mrqcof(atry,covar,da);
            if (fabs(chisq-ochisq) < std::max(tol,tol*chisq)) done++;
			if (chisq < ochisq) {
				alamda *= 0.1;
				ochisq=chisq;
				for (j=0;j<mfit;j++) {
					for (k=0;k<mfit;k++) alpha[j][k]=covar[j][k];
						beta[j]=da[j];
				}
				for (l=0;l<ma;l++) a[l]=atry[l];
			} else {
				alamda *= 10.0;
				chisq=ochisq;
			}
        }
        return 2; // too many iterations
	}

	void mrqcof(VecDoub_I &a, MatDoub_O &alpha, VecDoub_O &beta) {
        int i,j,k,l,m;
        double ymod,wt,sig2i,dy;
		VecDoub dyda(ma);
		for (j=0;j<mfit;j++) {
			for (k=0;k<=j;k++) alpha[j][k]=0.0;
			beta[j]=0.;
		}
		chisq=0.;
		for (i=0;i<ndat;i++) {
            funcs(x[i], x2[i],a,ymod,dyda);
			sig2i=1.0/(sig[i]*sig[i]);
			dy=y[i]-ymod;
			for (j=0,l=0;l<ma;l++) {
				if (ia[l]) {
					wt=dyda[l]*sig2i;
					for (k=0,m=0;m<l+1;m++)
						if (ia[m]) alpha[j][k++] += wt*dyda[m];
					beta[j++] += dy*wt;
				}
			}
			chisq += dy*dy*sig2i;
		}
		for (j=1;j<mfit;j++)
			for (k=0;k<j;k++) alpha[k][j]=alpha[j][k];
	}

	void covsrt(MatDoub_IO &covar) {
        int i,j,k;
		for (i=mfit;i<ma;i++)
			for (j=0;j<i+1;j++) covar[i][j]=covar[j][i]=0.0;
		k=mfit-1;
		for (j=ma-1;j>=0;j--) {
			if (ia[j]) {
                for (i=0;i<ma;i++) std::swap(covar[i][k],covar[i][j]);
                for (i=0;i<ma;i++) std::swap(covar[k][i],covar[j][i]);
				k--;
			}
		}
	}

    std::vector<std::vector<double>> plot(){
        std::vector<std::vector<double>> data{std::vector<double>{},std::vector<double>{}};
        double volume_step{},x1_stock{};
        volume_step=(volume.back()-volume[0])/999;
        x1_stock=(x.back()*volume.back()-x[0]*volume[0])/(volume.back()-volume[0]);
        double y_calc;
        for (int i=0;i<1000;i++){
            double x1{(x[0]*volume[0]+i*x1_stock*volume_step) / (volume[0]+i*volume_step)};
            double x_prot{(x2[0]*volume[0])/(volume[0]+i*volume_step)};
            y_calc=funcs2(x1,x_prot,a);
            data[1].push_back(y_calc);
            data[0].push_back(x1);
        }
        return data;
    }
};

struct Fitmrq2 {
    static const int NDONE=4, ITMAX=1000;
    int ndat, ma, mfit;
    VecDoub_I x, x2,x3, y, volume, sig;
    double tol;
    void (*funcs)(const double, const double, const double, VecDoub_I &, double &, VecDoub_O &);
    double (*funcs2)(const double, const double, const double,const std::vector<double> &);
    VecBool ia;
    MatDoub alpha;
    VecDoub a;
    MatDoub covar;
    double chisq;
    bool did_gauss_work;

    Fitmrq2(VecDoub_I &xx, VecDoub_I &xx2,VecDoub_I &xx3,VecDoub_I &yy, VecDoub_I & v,VecDoub_I &ssig, VecDoub_I &aa,
    void funks(const double, const double, const double,VecDoub_I &, double &, VecDoub_O &), double funks2(const double, const double, const double,const std::vector<double> &),
           const double TOL=1.e-3) : ndat(xx.size()), ma(aa.size()), x(xx), x2(xx2),x3(xx3),y(yy), volume(v), sig(ssig),
    tol(TOL), funcs(funks), funcs2(funks2),ia(ma), alpha(ma,ma), a(aa), covar(ma,ma) {
        for (int i=0;i<ma;i++)
            ia[i] = true;
    }

    Fitmrq2(VecDoub_I &xx, VecDoub_I &xx2,VecDoub_I &xx3,VecDoub_I &yy, VecDoub_I & v,VecDoub_I &ssig, VecDoub_I &aa, VecBool &iia,
    void funks(const double, const double, const double,VecDoub_I &, double &, VecDoub_O &), double funks2(const double, const double, const double,const std::vector<double> &),
           const double TOL=1.e-3) : ndat(xx.size()), ma(aa.size()), x(xx), x2(xx2),x3(xx3),y(yy), volume(v), sig(ssig),
    tol(TOL), funcs(funks), funcs2(funks2),ia(ma), alpha(ma,ma), a(aa), covar(ma,ma) {
    if (iia.size() == ia.size())
        ia = iia;
    else
        for (int i=0;i<ma;i++)
            ia[i] = true;
    }

    int fit() {
        int j,k,l,iter,done=0;
        double alamda=.001,ochisq;
        VecDoub atry(ma),beta(ma),da(ma);
        mfit=0;
        for (j=0;j<ma;j++) if (ia[j]) mfit++;
        MatDoub oneda(mfit,1), temp(mfit,mfit);
        mrqcof(a,alpha,beta);
        for (j=0;j<ma;j++) atry[j]=a[j];
        ochisq=chisq;
        for (iter=0;iter<ITMAX;iter++) {
            if (done==NDONE) alamda=0.;
            for (j=0;j<mfit;j++) {
                for (k=0;k<mfit;k++) covar[j][k]=alpha[j][k];
                covar[j][j]=alpha[j][j]*(1.0+alamda);
                for (k=0;k<mfit;k++) temp[j][k]=covar[j][k];
                oneda[j][0]=beta[j];
            }
            did_gauss_work=gaussj(temp,oneda);
            if (did_gauss_work==false)
                return 1;
            for (j=0;j<mfit;j++) {
                for (k=0;k<mfit;k++) covar[j][k]=temp[j][k];
                da[j]=oneda[j][0];
            }
            if (done==NDONE) {
                covsrt(covar);
                covsrt(alpha);
                return 0;
            }
            for (j=0,l=0;l<ma;l++)
                if (ia[l]) atry[l]=a[l]+da[j++];
            mrqcof(atry,covar,da);
            if (fabs(chisq-ochisq) < std::max(tol,tol*chisq)) done++;
            if (chisq < ochisq) {
                alamda *= 0.1;
                ochisq=chisq;
                for (j=0;j<mfit;j++) {
                    for (k=0;k<mfit;k++) alpha[j][k]=covar[j][k];
                        beta[j]=da[j];
                }
                for (l=0;l<ma;l++) a[l]=atry[l];
            } else {
                alamda *= 10.0;
                chisq=ochisq;
            }
        }
        return 2; // too many iterations
    }


    void mrqcof(VecDoub_I &a, MatDoub_O &alpha, VecDoub_O &beta) {
        int i,j,k,l,m;
        double ymod,wt,sig2i,dy;
        VecDoub dyda(ma);
        for (j=0;j<mfit;j++) {
            for (k=0;k<=j;k++) alpha[j][k]=0.0;
            beta[j]=0.;
        }
        chisq=0.;
        for (i=0;i<ndat;i++) {
            funcs(x[i], x2[i], x3[i],a,ymod,dyda);
            sig2i=1.0/(sig[i]*sig[i]);
            dy=y[i]-ymod;
            for (j=0,l=0;l<ma;l++) {
                if (ia[l]) {
                    wt=dyda[l]*sig2i;
                    for (k=0,m=0;m<l+1;m++)
                        if (ia[m]) alpha[j][k++] += wt*dyda[m];
                    beta[j++] += dy*wt;
                }
            }
            chisq += dy*dy*sig2i;
        }
        for (j=1;j<mfit;j++)
            for (k=0;k<j;k++) alpha[k][j]=alpha[j][k];
    }

    void covsrt(MatDoub_IO &covar) {
        int i,j,k;
        for (i=mfit;i<ma;i++)
            for (j=0;j<i+1;j++) covar[i][j]=covar[j][i]=0.0;
        k=mfit-1;
        for (j=ma-1;j>=0;j--) {
            if (ia[j]) {
                for (i=0;i<ma;i++) std::swap(covar[i][k],covar[i][j]);
                for (i=0;i<ma;i++) std::swap(covar[k][i],covar[j][i]);
                k--;
            }
        }
    }

    std::vector<std::vector<double>> plot(){
        std::vector<std::vector<double>> data{std::vector<double>{},std::vector<double>{}};
        double volume_step{},x1_stock{};
        volume_step=(volume.back()-volume[0])/999;
        x1_stock=(x.back()*volume.back()-x[0]*volume[0])/(volume.back()-volume[0]);
        double y_calc;
        for (int i=0;i<1000;i++){
            double x1{(x[0]*volume[0]+i*x1_stock*volume_step) / (volume[0]+i*volume_step)};
            double x_prot{(x2[0]*volume[0])/(volume[0]+i*volume_step)};
            double x_comp{x3[0]*volume[0]/(volume[0]+i*volume_step)};
            y_calc=funcs2(x1,x_prot,x_comp,a);
            data[1].push_back(y_calc);
            data[0].push_back(x1);
        }
        return data;
    }
};

struct Fitmrq3 {
    static const int NDONE=4, ITMAX=1000;
    int ndat, ma, mfit;
    VecDoub_I x, y, sig;
    double tol;
    void (*funcs)(const double, VecDoub_I &, double &, VecDoub_O &);
    double (*funcs2)(const double, const std::vector<double> &);
    VecBool ia;
    MatDoub alpha;
    VecDoub a;
    MatDoub covar;
    double chisq;
    bool did_gauss_work;

    Fitmrq3(VecDoub_I &xx, VecDoub_I &yy,VecDoub_I &ssig, VecDoub_I &aa,
    void funks(const double, VecDoub_I &, double &, VecDoub_O &), double funks2(const double, const std::vector<double> &),
           const double TOL=1.e-3) : ndat(xx.size()), ma(aa.size()), x(xx), y(yy), sig(ssig),
    tol(TOL), funcs(funks), funcs2(funks2),ia(ma), alpha(ma,ma), a(aa), covar(ma,ma) {
        for (int i=0;i<ma;i++)
            ia[i] = true;
    }

    Fitmrq3(VecDoub_I &xx, VecDoub_I &yy,VecDoub_I &ssig, VecDoub_I &aa, VecBool &iia,
    void funks(const double, VecDoub_I &, double &, VecDoub_O &), double funks2(const double, const std::vector<double> &),
           const double TOL=1.e-3) : ndat(xx.size()), ma(aa.size()), x(xx), y(yy), sig(ssig),
    tol(TOL), funcs(funks), funcs2(funks2),ia(ma), alpha(ma,ma), a(aa), covar(ma,ma) {
        if (iia.size() == ia.size())
            ia = iia;
        else
            for (int i=0;i<ma;i++)
                ia[i] = true;
    }

    void hold(const int i, const double val) {ia[i]=false; a[i]=val;}
    void free(const int i) {ia[i]=true;}

    int fit() {
        int j,k,l,iter,done=0;
        double alamda=.001,ochisq;
        VecDoub atry(ma),beta(ma),da(ma);
        mfit=0;
        for (j=0;j<ma;j++) if (ia[j]) mfit++;
        MatDoub oneda(mfit,1), temp(mfit,mfit);
        mrqcof(a,alpha,beta);
        for (j=0;j<ma;j++) atry[j]=a[j];
        ochisq=chisq;
        for (iter=0;iter<ITMAX;iter++) {
            if (done==NDONE) alamda=0.;
            for (j=0;j<mfit;j++) {
                for (k=0;k<mfit;k++) covar[j][k]=alpha[j][k];
                covar[j][j]=alpha[j][j]*(1.0+alamda);
                for (k=0;k<mfit;k++) temp[j][k]=covar[j][k];
                oneda[j][0]=beta[j];
            }
            did_gauss_work=gaussj(temp,oneda);
            if (did_gauss_work==false)
                return 1;
            for (j=0;j<mfit;j++) {
                for (k=0;k<mfit;k++) covar[j][k]=temp[j][k];
                da[j]=oneda[j][0];
            }
            if (done==NDONE) {
                covsrt(covar);
                covsrt(alpha);
                return 0;
            }
            for (j=0,l=0;l<ma;l++)
                if (ia[l]) atry[l]=a[l]+da[j++];
            mrqcof(atry,covar,da);
            if (fabs(chisq-ochisq) < std::max(tol,tol*chisq)) done++;
            if (chisq < ochisq) {
                alamda *= 0.1;
                ochisq=chisq;
                for (j=0;j<mfit;j++) {
                    for (k=0;k<mfit;k++) alpha[j][k]=covar[j][k];
                        beta[j]=da[j];
                }
                for (l=0;l<ma;l++) a[l]=atry[l];
            } else {
                alamda *= 10.0;
                chisq=ochisq;
            }
        }
        return 2;
    }

    void mrqcof(VecDoub_I &a, MatDoub_O &alpha, VecDoub_O &beta) {
        int i,j,k,l,m;
        double ymod,wt,sig2i,dy;
        VecDoub dyda(ma);
        for (j=0;j<mfit;j++) {
            for (k=0;k<=j;k++) alpha[j][k]=0.0;
            beta[j]=0.;
        }
        chisq=0.;
        for (i=0;i<ndat;i++) {
            funcs(x[i],a,ymod,dyda);
            sig2i=1.0/(sig[i]*sig[i]);
            dy=y[i]-ymod;
            for (j=0,l=0;l<ma;l++) {
                if (ia[l]) {
                    wt=dyda[l]*sig2i;
                    for (k=0,m=0;m<l+1;m++)
                        if (ia[m]) alpha[j][k++] += wt*dyda[m];
                    beta[j++] += dy*wt;
                }
            }
            chisq += dy*dy*sig2i;
        }
        for (j=1;j<mfit;j++)
            for (k=0;k<j;k++) alpha[k][j]=alpha[j][k];
    }

    void covsrt(MatDoub_IO &covar) {
        unsigned int i,j,k;
        for (i=unsigned(mfit);i<unsigned(ma);i++)
            for (j=0;j<i+1;j++) covar[i][j]=covar[j][i]=0.0;
        k=mfit-1;
        for (j=ma-1;j>=0;j--) {
            if (ia[j]) {
                for (i=0;i<ma;i++) std::swap(covar[i][k],covar[i][j]);
                for (i=0;i<ma;i++) std::swap(covar[k][i],covar[j][i]);
                k--;
            }
        }
    }

    std::vector<std::vector<double>> plot(){
        std::vector<std::vector<double>> data{std::vector<double>{},std::vector<double>{}};
        double y_calc;
        for (unsigned int i=0;i<x.size();i++){
            y_calc=funcs2(x[i], a);
            data[1].push_back(y_calc);
            data[0].push_back(x[i]);
        }
        return data;
    }
};
#endif
