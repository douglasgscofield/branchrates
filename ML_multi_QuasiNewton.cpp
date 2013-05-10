#include "ML_multi_QuasiNewton.h"


void ML_multi_QuasiNewton::dfpmin(v_ratep_type& p,
                                  const double gtol,
                                  int& iter,
                                  double& fret,
                                  ptr_eval_func func,
                                  ptr_eval_gradient_func dfunc
                                  // double func(const v_ratep_type& p),
                                  // void dfunc(const v_ratep_type& in,
                                  //            v_ratep_type& out)
                                 )
{
    // from Press et al 2002
    const int ITMAX = 200;
    const double EPS=std::numeric_limits<double>::epsilon();
    const double TOLX = 4.0*EPS, STPMX = 100.0;
    bool check;
    int i, its, j;
    double den, fac, fad, fae, fp, stpmax, sum=0.0, sumdg, sumxi, temp, test;

    int n = p.size();
    v_ratep_type dg(n), g(n), hdg(n), pnew(n), xi(n);
    std::vector<v_ratep_type> hessin;
    hessin.resize(n);
    for (int i = 0; i < n; ++i) hessin[i].resize(n);
    fp = (this->*func)(p);  // fp = func(p);
    (this->*dfunc)(p, g);   // dfunc(p, g);
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) hessin[i][j] = 0.0;
        hessin[i][i] = 1.0;
        xi[i] = -g[i];
        sum += p[i]*p[i];
    }
    stpmax = STPMX * MAX(std::sqrt(sum), double(n));
    for (its = 0; its < ITMAX; ++its) {
        iter = its;
        lnsrch(p, fp, g, xi, pnew, fret, stpmax, check, func);
        fp = fret;
        for (i = 0; i < n; ++i) {
            xi[i] = pnew[i] - p[i];
            p[i] = pnew[i];
        }
        test = 0.0;
        for (i = 0; i < n; ++i) {
            temp = std::abs(xi[i]) / MAX(std::abs(p[i]), 1.0);
            if (temp > test) test = temp;
        }
        if (test < TOLX)
            return;
        for (i = 0; i < n; ++i) dg[i] = g[i];
        (this->*dfunc)(p, g);  // dfunc(p, g);
        test = 0.0;
        den = MAX(fret, 1.0);
        for (i = 0; i < n; ++i) {
            temp = std::abs(g[i])*MAX(std::abs(p[i]), 1.0)/den;
            if (temp > test) test = temp;
        }
        if (test < gtol)
            return;
        for (i = 0; i < n; ++i) dg[i] = g[i] - dg[i];
        for (i = 0; i < n; ++i) {
            hdg[i] = 0.0;
            for (j = 0; j < n; ++j) hdg[i] += hessin[i][j]*dg[j];
        }
        fac = fae = sumdg = sumxi = 0.0;
        for (i = 0; i < n; ++i) {
            fac += dg[i]*xi[i];
            fae += dg[i]*hdg[i];
            sumdg += SQR(dg[i]);
            sumxi += SQR(xi[i]);
        }
        if (fac > std::sqrt(EPS*sumdg*sumxi)) {
            fac = 1.0 / fac;
            fad = 1.0 / fae;
            for (i = 0; i < n; ++i) dg[i] = fac*xi[i] - fad*hdg[i];
            for (i = 0; i < n; ++i) {
                for (j = i; j < n; ++j) {
                    hessin[i][j] += fac*xi[i]*xi[j] -
                                    fad*hdg[i]*hdg[j] +
                                    fae*dg[i]*dg[j];
                    hessin[j][i] = hessin[i][j];
                }
            }
        }
        for (i = 0; i < n; ++i) {
            xi[i] = 0.0;
            for (j = 0; j < n; ++j) xi[i] -= hessin[i][j]*g[j];
        }
    }
    std::cerr << "dfpmin(): too many iterations" << std::endl;
    assert(false);
}


