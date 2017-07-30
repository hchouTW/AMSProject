#ifndef __CPPLibs_MGMath_C__
#define __CPPLibs_MGMath_C__


namespace MGMath {


inline double Quantile(double prob) {
    // Computes quantiles for standard normal distribution N(0, 1)
    // at probability prob
    // ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, 477-484.

    if ((prob<=0)||(prob>=1)) {
        std::cerr << "MGMath::Quantile : probability outside (0, 1)\n";
        return 0;
    }

    constexpr double a0 = 3.3871328727963666080e0;
    constexpr double a1 = 1.3314166789178437745e+2;
    constexpr double a2 = 1.9715909503065514427e+3;
    constexpr double a3 = 1.3731693765509461125e+4;
    constexpr double a4 = 4.5921953931549871457e+4;
    constexpr double a5 = 6.7265770927008700853e+4;
    constexpr double a6 = 3.3430575583588128105e+4;
    constexpr double a7 = 2.5090809287301226727e+3;
    constexpr double b1 = 4.2313330701600911252e+1;
    constexpr double b2 = 6.8718700749205790830e+2;
    constexpr double b3 = 5.3941960214247511077e+3;
    constexpr double b4 = 2.1213794301586595867e+4;
    constexpr double b5 = 3.9307895800092710610e+4;
    constexpr double b6 = 2.8729085735721942674e+4;
    constexpr double b7 = 5.2264952788528545610e+3;
    constexpr double c0 = 1.42343711074968357734e0;
    constexpr double c1 = 4.63033784615654529590e0;
    constexpr double c2 = 5.76949722146069140550e0;
    constexpr double c3 = 3.64784832476320460504e0;
    constexpr double c4 = 1.27045825245236838258e0;
    constexpr double c5 = 2.41780725177450611770e-1;
    constexpr double c6 = 2.27238449892691845833e-2;
    constexpr double c7 = 7.74545014278341407640e-4;
    constexpr double d1 = 2.05319162663775882187e0;
    constexpr double d2 = 1.67638483018380384940e0;
    constexpr double d3 = 6.89767334985100004550e-1;
    constexpr double d4 = 1.48103976427480074590e-1;
    constexpr double d5 = 1.51986665636164571966e-2;
    constexpr double d6 = 5.47593808499534494600e-4;
    constexpr double d7 = 1.05075007164441684324e-9;
    constexpr double e0 = 6.65790464350110377720e0;
    constexpr double e1 = 5.46378491116411436990e0;
    constexpr double e2 = 1.78482653991729133580e0;
    constexpr double e3 = 2.96560571828504891230e-1;
    constexpr double e4 = 2.65321895265761230930e-2;
    constexpr double e5 = 1.24266094738807843860e-3;
    constexpr double e6 = 2.71155556874348757815e-5;
    constexpr double e7 = 2.01033439929228813265e-7;
    constexpr double f1 = 5.99832206555887937690e-1;
    constexpr double f2 = 1.36929880922735805310e-1;
    constexpr double f3 = 1.48753612908506148525e-2;
    constexpr double f4 = 7.86869131145613259100e-4;
    constexpr double f5 = 1.84631831751005468180e-5;
    constexpr double f6 = 1.42151175831644588870e-7;
    constexpr double f7 = 2.04426310338993978564e-15;

    constexpr double split1 = 0.425;
    constexpr double split2=5.;
    constexpr double konst1=0.180625;
    constexpr double konst2=1.6;

    double q, r, quantile;
    q=prob-0.5;
    if (std::fabs(q)<split1) {
        r=konst1-q*q;
        quantile = q* (((((((a7 * r + a6) * r + a5) * r + a4) * r + a3)
                        * r + a2) * r + a1) * r + a0) /
            (((((((b7 * r + b6) * r + b5) * r + b4) * r + b3)
               * r + b2) * r + b1) * r + 1.);
    } else {
        if(q<0) r=prob;
        else    r=1-prob;
        //error case
        if (r<=0)
            quantile=0;
        else {
            r=std::sqrt(-std::log(r));
            if (r<=split2) {
                r=r-konst2;
                quantile=(((((((c7 * r + c6) * r + c5) * r + c4) * r + c3)
                                * r + c2) * r + c1) * r + c0) /
                    (((((((d7 * r + d6) * r + d5) * r + d4) * r + d3)
                       * r + d2) * r + d1) * r + 1);
            } else{
                r=r-split2;
                quantile=(((((((e7 * r + e6) * r + e5) * r + e4) * r + e3)
                                * r + e2) * r + e1) * r + e0) /
                    (((((((f7 * r + f6) * r + f5) * r + f4) * r + f3)
                       * r + f2) * r + f1) * r + 1);
            }
            if (q<0) quantile=-quantile;
        }
    }
    return quantile;
}


inline double Quantile(double prob, double wgtA, double sgmA, double wgtB, double sgmB) {
    double _wgtA = wgtA / (wgtA + wgtB);
    double _wgtB = wgtB / (wgtA + wgtB);
    double _sgmA = (sgmA / sgmA);
    double _sgmB = (sgmB / sgmA);

    return 0.0;
}


}


#endif // __CPPLibs_MGMath_C__
