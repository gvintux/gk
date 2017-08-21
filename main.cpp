#include <iostream>
#include <iomanip>
#include <cfloat>

using namespace std;
typedef double (*integrand_1d) (double, void*);

const double gk15_roots [] =
{
    0,
    0.20778495500789846760068940377324491348,
    0.40584515137739716690660641207696146335,
    0.58608723546769113029414483825872959844,
    0.74153118559939443986386477328078840707,
    0.8648644233597690727897127886409262012,
    0.9491079123427585245261896840478512624,
    0.99145537112081263920685469752632851664
};

const double gk15_weights [] =
{
    0.2094821410847278280129991748917142637,
    0.20443294007529889241416199923464908472,
    0.1903505780647854099132564024210136828,
    0.1690047266392679028265834265985502841,
    0.1406532597155259187451895905102379204,
    0.1047900103222501838398763225415180174,
    0.0630920926299785532907006631892042867,
    0.022935322010529224963732008058969592
};
const double g7_weights [] =
{
    0.41795918367346938775510204081632653061,
    0.38183005050511894495036977548897513388,
    0.2797053914892766679014677714237795825,
    0.1294849661688696932706114326790820183
};

double fun(double x, void *args)
{
    return x*x;
}

void integrate_1d_gk(const int n, const double gk_roots[], const double g_weights[], const double gk_weights[],
                     const integrand_1d f, void *args, double a, double b, double *result, double *abserr,
                     double *resabs, double *resasc)
{
    const double half = (b-a) / 2.0;
    const double middle = (a + b) / 2.0;
    double k_result = gk_weights[0] * f(gk_roots[0] * half + middle, args);
    for(int i = 1; i <= n / 2; ++i) k_result += gk_weights[i] * (f(gk_roots[i] * half + middle, args) +  f(-gk_roots[i] * half + middle, args));
    k_result *= half;
    (*result) = k_result;
}

int main()
{
    double res;
    integrate_1d_gk(15, gk15_roots, g7_weights, gk15_weights, fun, 0, 1, 3, &res, 0, 0, 0);
    cout << setprecision(DBL_DIG);
    cout << res << endl;
    return 0;
}