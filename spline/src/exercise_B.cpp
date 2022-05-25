#include"Spline.hpp"

using namespace std;

void B()
{
    cout << "================ Assignment.B ===============" << endl;
    function<double(double)> func[2]={[](double x){return 1.0/(1+25*x*x);},
                                      [](double x){return -50*x/pow((25*x*x+1),2);}};
    int N[5] = {5,10,20,40,80};
    for(auto n : N)
    {
        double a[n+1];
        int b[n+1];
        for(int i=0; i<n+1; i++){
            a[i]=-1+(2.0/n)*i;
            if(i>0 && i<n) b[i]=1;
            else b[i]=2;
        }
        InterpConditions inter(n+1,a,b,func);
        auto res1 = interpolate<4>(inter,complete);
        auto res2 = interpolate<4>(inter,notAknot);
        cout << "complete cubic spline:" << endl << res1;
        cout << "not-A-knot cubic spline:" << endl << res2;
        cout << endl;
        cout << "When n = " << n << ":" << endl;
        cout << "Error for complete cubic spline:" << endl;
        double error = 0;
        for(int i=0; i<n; i++){
            double dx = 0.5*(a[i]+a[i+1]);
            double val = res1(dx)-func[0](dx);
            if(abs(val) > error)
                error = val;
        }
        cout << "   " << error << endl;
        cout << "Error for notAknot cubic spline:" << endl;
        error = 0;
        for(int i=0; i<n; i++){
            double dx = 0.5*(a[i]+a[i+1]);
            double val = res2(dx)-func[0](dx);
            if(abs(val) > error)
                error = val;
        }
        cout << "   " << error << endl << endl;
    }

}

int main(int argc, char* argv[])
{
    B();
    return 0;
}