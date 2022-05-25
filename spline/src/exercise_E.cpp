#include"Spline.hpp"
#include<iomanip>

using namespace std;

void E()
{
    cout << "================= Assignment.E =================" << endl;
    int list[3] = {10,40,160};
    cout.precision(12);
#define pi M_PI
    for(auto n : list){
        double t[n+1];
        double r[n+1];
        vector<Vec<double,2>> vec(n+1);
        for(int i=0; i<n+1; i++){
            t[i] = 2*pi*i*1.0/n;
            r[i] = sqrt( 3/( sin(t[i])*( 0.25*sin(t[i])-3*abs(cos(t[i])) ) +2) );
            vec[i] = (Vec<double,2>){r[i]*cos(t[i]),r[i]*sin(t[i])};
        }
        auto res = fitCurve<4>(vec,periodic);
        cout << "when n = " << n << ":" << endl;
        output(res);
    }
#undef pi
}

int main(int argc, char* argv[])
{
    E();
    return 0;
}