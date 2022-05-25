#include"Spline.hpp"

using namespace std;

void F()
{
    cout << "================= Assignment.F =================" << endl;
    cout.precision(12);
    double a1[9]={1,2,5,6,7,8,10,13,17};
    int b1[9]={2,1,1,1,1,1,1,1,2};
    double c1[11]={3,1,3.7,3.9,4.2,5.7,6.6,7.1,6.7,4.5,-0.67};
    InterpConditions inter1(9,a1,b1,c1);
    auto res1 = interpolate<4>(inter1,complete);
    cout << "First piece:" << endl << res1 << endl;
    double a2[7]={17,20,23,24,25,27,27.7};
    int b2[7]={2,1,1,1,1,1,2};
    double c2[9]={4.5,3,7,6.1,5.6,5.8,5.2,4.1,-4};
    InterpConditions inter2(7,a2,b2,c2);
    auto res2 = interpolate<4>(inter2,complete);
    cout << "Second piece: " << endl << res2 << endl;
    double a3[4]={27.7,28,29,30};
    int b3[4]={2,1,1,2};
    double c3[6]={4.1,0.33,4.3,4.1,3,-1.5};
    InterpConditions inter3(4,a3,b3,c3);
    auto res3 = interpolate<4>(inter3,complete);
    cout << "Third piece:" << endl << res3;

}

int main(int argc, char* argv[])
{
    F();
    return 0;
}