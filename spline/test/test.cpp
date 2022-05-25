#include"Spline.hpp"

using namespace std;

int main()
{
//=============== Test for Polynomial =================
    //测试基本的多项式声明
    Polynomial<3,int> poly{1,2,3,4};
    cout << poly(1) << endl;
    Polynomial<3,double> p2;
    p2 = poly;
    Polynomial<2,double> p3{5,7,9};

    //测试多项式加减法
    auto pp = p2 + poly;
    auto pp2 = p2 - p3;
    cout << pp[0] << " " << pp[1] << " " << pp[2] << endl;
    cout << pp2[0] << " " << pp2[1] << " " << pp2[2] << endl;

    //测试“加一个数”的重载方法
    auto px = p2 + 5;
    cout << px << endl;

    //测试多项式乘法
    auto p4 = p2*p3;
    cout << p4 << endl;

    //测试求导
    auto p5 = p2.diff();
    cout << p5 << endl;

    //测试系数为向量的多项式及其基本功能
    Polynomial<3,Vec<double,2>> pol{{0,1},{0,1},{1,2},{1,2}};
    cout << pol << endl;
    Polynomial<2,Vec<double,2>> pol2{{2,3},{5,5},{1,0}};
    cout << pol+pol2 << endl;
    cout << pol*pol2 << endl;
    cout << pol(1) << endl;

//============ Test for one-dimension ppForm spline =============
    function<double(double)> func[2] = {[](double x){return 1/x;},
                                        [](double x){return -1/(x*x);}};
    double a[4] = {1,2,3,4};
    int b[4] = {2,1,1,2};
    InterpConditions inter(4,a,b,func);

    //分别测试三种不同的ppForm样条
    auto RES = interpolate<4>(inter,complete);
    auto RES2 = interpolate<4>(inter,notAknot);
    auto RES3 = interpolate<4>(inter,periodic);
    cout << RES << RES2 << RES3 ;
    //验证生成的样条是对的
    cout << RES3.spline[0].diff().diff()(1) << " " << RES3.spline[2].diff().diff()(4) << endl;
    auto RES4 = interpolate<2>(inter,complete);
    cout << RES4;

    

//========== Test for two-dimension ppForm spline ==========
    //测试fitCurve方法对order=4的complete情形
    vector<Vec<double,2>> lst{{1,1},{4,5},{10,13},{1,1},{1,1}};
    auto RES5 = fitCurve<4>(lst,complete);
    cout << RES5;
    cout << RES5.spline[0](5) <<" " << RES5.spline[1](5)<< endl;
    
    //测试fitCurve方法对order=4的notAknot情形
    vector<Vec<double,2>> lst2{{1,1},{7,9},{4,5},{9,17}};
    auto RES6 = fitCurve<4>(lst2,notAknot);
    cout << RES6;
    
    //测试fitCurve方法对order=4的periodic情形
    auto RES7 = fitCurve<4>(lst2,periodic);
    cout << RES7;
    // cout << RES7.spline[0].diff()(0) << " " << RES7.spline[2].diff()(28) << endl;
    // cout << RES7.spline[0].diff().diff()(0) << " " << RES7.spline[2].diff().diff()(28) << endl;
    
    //测试fitCurve方法对order=2的情形
    auto RES8 = fitCurve<2>(lst2,complete);
    cout << RES8;


//================ Test for cardinalB spline. ================
    function<double(double)> _func[2] = {[](double x){return 1/(1+x*x);},
                                         [](double x){return -2*x/((1+x*x)*(1+x*x));}};
    
    //测试Corollary 4.58
    double _a[11] = {-5,-4,-3,-2,-1,0,1,2,3,4,5};
    int _b[11] = {2,1,1,1,1,1,1,1,1,1,2};
    InterpConditions _interp(11,_a,_b,_func);
    auto RRRES = interpolate<3>(_interp);
    cout << RRRES;

    //测试Corollary 4.59
    #define f _func[0]
    double _a_[10] = {-5,-4,-3,-2,-1,0,1,2,3,4};
    int _b_[10] = {2,1,1,1,1,1,1,1,1,2};
    double _c[12] = {f(-4.5),f(-5),f(-3.5),f(-2.5),f(-1.5),
                    f(-0.5),f(0.5),f(1.5),f(2.5),f(3.5),f(4.5),f(5)};
    InterpConditions _interp2(10,_a_,_b_,_c);
    auto RRES2 = interpolate<2>(_interp2);
    cout << RRES2;
    #undef f

    return 0;
}