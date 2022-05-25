#include"Spline.hpp"

using namespace std;

void C_D()
{
    cout << "============== Assignment.C ==============" << endl;
    cout.precision(12);
    function<double(double)> _func[2] = {[](double x){return 1/(1+x*x);},
                                         [](double x){return -2*x/((1+x*x)*(1+x*x));}};
    double _a[11] = {-5,-4,-3,-2,-1,0,1,2,3,4,5};
    int _b[11] = {2,1,1,1,1,1,1,1,1,1,2};
    InterpConditions _interp(11,_a,_b,_func);
    auto RES = interpolate<3>(_interp);
    cout << "Result for Corollary 4.58:" << endl;
    cout << RES << endl;
#define f _func[0]
    double _a_[10] = {-5,-4,-3,-2,-1,0,1,2,3,4};
    int _b_[10] = {2,1,1,1,1,1,1,1,1,2};
    double _c[12] = {f(-4.5),f(-5),f(-3.5),f(-2.5),f(-1.5),
                    f(-0.5),f(0.5),f(1.5),f(2.5),f(3.5),f(4.5),f(5)};
    InterpConditions _interp2(10,_a_,_b_,_c);
    auto RES2 = interpolate<2>(_interp2);
    cout << "Result for Corollary 4.59:" << endl;
    cout << RES2 << endl;

    cout << "============== Assignment.D ==============" << endl;
    double x[7] = {-3.5,-3,-0.5,0,0.5,3,3.5};
    cout << "For the corollary 4.58:" << endl;
    for(auto &i : x)
        cout << "   E(" << i << ") = " << abs(RES(i)-f(i)) << endl;
    cout << "For the corollary 4.59:" << endl;
    for(auto &i : x)
        cout << "   E(" << i << ") = " << abs(RES2(i)-f(i)) << endl;
#undef f
}

int main(int argc, char* argv[])
{
    C_D();
    return 0;
}