#include "ML_multi_Powell.h"


#define   DEBUG_LINMIN   0

void ML_multi_Powell::powell(v_ratep_type& p,
                             mat_ratep_type& xi,
                             const double ftol,
                             int& iter,
                             double& fret,
                             ptr_eval_func func
                             )
{
    const int ITMAX = 200;
    const double TINY = 1.0e-25;
    int i, j, ibig;
    double del, fp, fptt, t;

    int n = p.size();
    v_ratep_type pt(n), ptt(n), xit(n);
    fret = (this->*func)(p);
    for (j = 0; j < n; ++j) pt[j] = p[j];
    for (iter = 0; ; ++iter) {
        fp = fret;
        ibig = 0;
        del = 0.0;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j) xit[j] = xi[j][i];
            fptt = fret;
            linmin(p, xit, fret, func);
            if (fptt - fret > del) {
                del = fptt - fret;
                ibig = i + 1;
            }
        }
        if (2.0*(fp - fret) <= ftol*(std::abs(fp) + std::abs(fret)) + TINY) {
            return;
        }
        if (iter == ITMAX) {
            std::cerr << "powell(): exceeded max iterations " << ITMAX << std::endl;
            assert(false);
        }
        for (j = 0; j < n; ++j) {
            ptt[j] = 2.0*p[j] - pt[j];
            xit[j] = p[j] - pt[j];
            pt[j] = p[j];
        }
        fptt = (this->*func)(ptt);
        if (fptt < fp) {
            t = 2.0*(fp - 2.0*fret + fptt)*SQR(fp - fret - del) - del*SQR(fp - fptt);
            if (t < 0.0) {
                linmin(p, xit, fret, func);
                for (j = 0; j < n; ++j) {
                    xi[j][ibig-1] = xi[j][n-1];
                    xi[j][n-1] = xit[j];
                }
            }
        }
    }
}


void ML_multi_Powell::linmin(v_ratep_type& p,
                             v_ratep_type& xi,
                             double& fret,
                             ptr_eval_func func
                             )
{
    int j;
    const double TOL=1.0e-8;
    double xx, xmin, fx, fb, fa, bx, ax;
    int n = p.size();
    ncom = n;
    pcom.resize(n);
    xicom.resize(n);
    nrfunc = func;
    for(j = 0; j < n; ++j) {
        pcom[j] = p[j];
        xicom[j] = xi[j];
    }
    ax = 0.0;
    xx = 1.0;
    mnbrak(ax, xx, bx, fa, fx, fb, &ML_multi_Powell::f1dim);
    fret = brent(ax, xx, bx, &ML_multi_Powell::f1dim, TOL, xmin);
    for (j = 0; j < n; ++j) {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
}


double ML_multi_Powell::f1dim(const double x)
{
    int j;
    v_ratep_type xt(ncom);
    for (j = 0; j < ncom; ++j) {
        xt[j] = pcom[j] + x*xicom[j];
    }
    return ((this->*nrfunc)(xt));
}


void ML_multi_Powell::mnbrak(double& ax,
                             double& bx,
                             double& cx,
                             double& fa,
                             double& fb,
                             double& fc,
                             ptr_eval_1_func f
                            )
{
    const double GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-20;
    double ulim, u, r, q, fu;
    fa = (this->*f)(ax);
    fb = (this->*f)(bx);
    if (fb > fa) {
        SWAP(ax, bx);
        SWAP(fb, fa);
    }
    cx = bx + GOLD*(bx - ax);
    fc = (this->*f)(cx);
    while (fb > fc) {
        r = (bx - ax)*(fb - fc);
        q = (bx - cx)*(fb - fa);
        u = bx - ((bx - cx)*q - (bx - ax)*r) / 
                 (2.0*SIGN(MAX(std::abs(q - r), TINY), q - r));
        ulim = bx + GLIMIT*(cx - bx);
        if ((bx - u)*(u - cx) > 0.0) {
            fu = (this->*f)(u);
            if (fu < fc) {
                ax = bx;
                bx = u;
                fa = fb;
                fb = fu;
                return;
            } else if (fu > fb) {
                cx = u;
                fc = fu;
                return;
            }
            u = cx + GOLD*(cx - bx);
            fu = (this->*f)(u);
        } else if ((cx - u)*(u - ulim) > 0.0) {
            fu = (this->*f)(u);
            if (fu < fc) {
                shft3(bx, cx, u, u + GOLD*(u - cx));
                shft3(fb, fc, fu, (this->*f)(u));
            }
        } else if ((u - ulim)*(ulim - cx) >= 0.0) {
            u = ulim;
            fu = (this->*f)(u);
        } else {
            u = cx + GOLD*(cx - bx);
            fu = (this->*f)(u);
        }
        shft3(ax, bx, cx, u);
        shft3(fa, fb, fc, fu);
    }
}


double ML_multi_Powell::brent(const double ax,
                              const double bx,
                              const double cx,
                              ptr_eval_1_func f,
                              const double tol,
                              double& xmin
                             )
{
    const int ITMAX = 100;
    const double CGOLD = 0.3819660;
    const double ZEPS = std::numeric_limits<double>::epsilon()*1.0e-3;

    int iter;
    double a, b, d = 0.0, etemp, fu, fv, fw, fx;
    double p, q, r, tol1, tol2, u, v, w, x, xm;
    double e = 0.0;

    a = (ax < cx ? ax : cx);
    b = (ax > cx ? ax : cx);
    x = w = v = bx;
    fw = fv = fx = (this->*f)(x);
    for (iter = 0; iter < ITMAX; ++iter) {
        xm = 0.5*(a + b);
        tol2 = 2.0*(tol1 = tol*std::abs(x) + ZEPS);
        if (std::abs(x - xm) <= (tol2 - 0.5*(b - a))) {
            xmin = x;
            return(fx);
        }
        if (std::abs(e) > tol1) {
            r = (x - w)*(fx - fv);
            q = (x - v)*(fx - fw);
            p = (x - v)*q - (x - w)*r;
            q = 2.0*(q - r);
            if (q > 0.0) p = -p;
            q = std::abs(q);
            etemp = e;
            e = d;
            if (std::abs(p) >= std::abs(0.5*q*etemp) || p <= q*(a - x) 
                || p >= q*(b - x)) {
                d = CGOLD*(e = (x >= xm ? a - x : b - x));
            } else {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = SIGN(tol1, xm - x);
            }
        } else {
            d = CGOLD*(e = (x >= xm ? a - x : b - x));
        }
        u = (std::abs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
        fu = (this->*f)(u);
        if (fu <= fx) {
            if (u >= x) a = x; else b = x;
            shft3(v, w, x, u);
            shft3(fv, fw, fx, fu);
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }
    }
    std::cerr << "brent(): exceeded iteration limit " << ITMAX << std::endl;
    assert(false);
    xmin = x;
    return fx;
}


