#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>

using namespace std;
typedef double (*integrand_1d) (double, void*);
typedef double (*integrand_2d) (double, double, void*);


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
    return 1/(1+x*x);
}

double fun2(double x, double y, void *args)
{
    return x+y;
}

void integrate_1d_gk(const int n, const double gk_roots[], const double g_weights[], const double gk_weights[],
                     const integrand_1d f, void *args, double a, double b, double *result, double *abserr,
                     double *resabs, double *resasc)
{
    const double half = (b-a) / 2.0;
    const double middle = (a + b) / 2.0;
    double f_middle = f(gk_roots[0] * half + middle, args);
    double k_result = gk_weights[0] * f_middle;
    double g_result = g_weights[0] * f_middle;

    for(int i = 1; i < n / 4 + 1; ++i)
    {
        double f_p = f(gk_roots[i * 2] * half + middle, args);
        double f_m = f(-gk_roots[i * 2] * half + middle, args);
        g_result += g_weights[i] * (f_p + f_m);
        k_result += gk_weights[i * 2] * (f_p + f_m);
    }

    for(int i = 0; i < n / 4 + 1; ++i)
    {
        double f_p = f(gk_roots[2 * i + 1] * half + middle, args);
        double f_m = f(-gk_roots[2 * i + 1] * half + middle, args);
        k_result += gk_weights[2 * i + 1] * (f_p + f_m);
    }
    k_result *= half;
    g_result *= half;
    (*result) = k_result;
    (*abserr) = pow(200 * fabs(g_result - k_result), 1.5);
}

void integrate_2d_gk(const int n, const double gk_roots[], const double g_weights[], const double gk_weights[],
                     const integrand_2d f, void *args, double a, double b, double c, double d, double *result, double *abserr,
                     double *resabs, double *resasc)
{
    const double half_x = (b-a) / 2.0;
    const double middle_x = (a + b) / 2.0;
    const double half_y = (d-c) / 2.0;
    const double middle_y = (c + d) / 2.0;
    double f_middle = f(gk_roots[0] * half_x + middle_x, gk_roots[0] * half_y + middle_y, args);
    double k_result = gk_weights[0] * gk_weights[0] * f_middle;
    double g_result = g_weights[0] * g_weights[0] * f_middle;
    for(int i = 1; i < n / 4 + 1; ++i)
    {
        double x = gk_roots[i * 2] * half_x + middle_x;
        for(int j = 1; j < n / 4 + 1; ++j)
        {
            double y = gk_roots[j * 2] * half_x + middle_x;
            double f_pp = f(x, y, args);
            double f_pm = f(x, -y, args);
            double f_mp = f(-x, y, args);
            double f_mm = f(-x, -y, args);
            g_result += g_weights[i] * g_weights[j] * (f_pp + f_pm + f_mp + f_mm);
            k_result += gk_weights[i * 2] * gk_weights[j * 2] * (f_pp + f_pm + f_mp + f_mm);
        }
    }

    for(int i = 0; i < n / 4 + 1; ++i)
    {
        double x = gk_roots[2*i + 1] * half_x + middle_x;
        for(int j = 0; j < n / 4 + 1; ++j)
        {
            double y = gk_roots[2*j + 1] * half_y + middle_y;
            double f_pp = f(x, y, args);
            double f_pm = f(x, -y, args);
            double f_mp = f(-x, y, args);
            double f_mm = f(-x, -y, args);
            k_result += gk_weights[2*i + 1] * gk_weights[2*j + 1] * (f_pp + f_pm + f_mp + f_mm);
        }
    }
    k_result *= half_x * half_y;
    g_result *= half_x * half_y;
    (*result) = g_result;
    (*abserr) = pow(200 * fabs(g_result - k_result), 1.5);
}

int main()
{
    double res, abserr;
    integrate_2d_gk(15, gk15_roots, g7_weights, gk15_weights, fun2, 0, 0, 1, 0, 1,  &res, &abserr, 0, 0);
    cout << setprecision(DBL_DIG*2);
    cout << res << endl;
    cout << abserr << endl;
    cout << (abserr > DBL_EPSILON) << endl;
    return 0;
}