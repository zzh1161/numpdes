/**************************************************************
 > @Description    : main.cpp
 > @Version        : 1.0
 > @Author         : zhang-zh
 > @Date           : 2022-03-29 12:43
 > @LastEditTime   : 2022-03-30 20:24
**************************************************************/
#include<fstream>
#include<string>
#include"../include/BVP.hpp"
#include"../include/json.hpp"

using json = nlohmann::json;

const char *boundarystream[] = {
    "./src/boundaryCondi/boundary1.json",
    "./src/boundaryCondi/boundary2.json",
    "./src/boundaryCondi/boundary3.json"
};

#define GET_THE_COEFFICIENT                                               \
    json j;                                                               \
    std::ifstream jfile(boundarystream[i]);                               \
    jfile >> j;                                                           \
    jfile.close();                                                        \
    double a1 = j.at("BoundaryOfSquare").at("left").at("CoeffiOfNeum");   \
    double b1 = j.at("BoundaryOfSquare").at("left").at("CoeffiOfDirch");  \
    double a2 = j.at("BoundaryOfSquare").at("right").at("CoeffiOfNeum");  \
    double b2 = j.at("BoundaryOfSquare").at("right").at("CoeffiOfDirch"); \
    double a3 = j.at("BoundaryOfSquare").at("up").at("CoeffiOfNeum");     \
    double b3 = j.at("BoundaryOfSquare").at("up").at("CoeffiOfDirch");    \
    double a4 = j.at("BoundaryOfSquare").at("down").at("CoeffiOfNeum");   \
    double b4 = j.at("BoundaryOfSquare").at("down").at("CoeffiOfDirch");  \
    double a5 = j.at("BoundaryOfDisk").at("CoeffiOfNeum");                \
    double b5 = j.at("BoundaryOfDisk").at("CoeffiOfDirch");

/**************************************************************
 > @description: The function is exp(y+sin(x)), 
 >               and the disk is (x-2/3)^2+(y-1/2)^2=(1/4)^2
**************************************************************/
void for_the_first_test()
{
    std::cout << "!! FOR THE 1ST TESTING FUNCTION: " << std::endl
              << "///////////////////////////////////////////////////"
              << std::endl << std::endl;

    auto left_mixed = [](double a, double b){
        return [a,b](double x, double y){
            return -cos(x)*exp(y+sin(x))*a + exp(y+sin(x))*b;
        };
    };
    auto right_mixed = [](double a, double b){
        return [a,b](double x, double y){
            return cos(x)*exp(y+sin(x))*a + exp(y+sin(x))*b;
        };
    };
    auto upper_mixed = [](double a, double b){
        return [a,b](double x, double y){
            return exp(y+sin(x))*(a+b);
        };
    };
    auto lower_mixed = [](double a, double b){
        return [a,b](double x, double y){
            return exp(y+sin(x))*(-a+b);
        };
    };
    auto circle_mixed = [](double a, double b){
        return [a,b](double x, double y){
            double x0 = 2.0/3, y0 = 0.5;
            double r = 0.25;
            return exp(y+sin(x))*(a*(x-x0)*cos(x)/r + a*(y-y0)/r + b);
        };
    };

    std::cout << "! For The Regular Domin: " << std::endl 
              << std::endl;
    for(int i=0; i<3; i++){
        GET_THE_COEFFICIENT
        BoundaryCondi<domin::square> B(left_mixed(a1,b1),
                                       right_mixed(a2,b2),
                                       upper_mixed(a3,b3),
                                       lower_mixed(a4,b4),
                                       j);
        std::vector<int> N = {8,16,32,64};
        std::cout << "No." << i+1 << std::endl;
        for(auto n: N){
            grid<domin::square> test(n);
            test.BVP_solver(B,[](double x, double y){return (sin(x)-1-cos(x)*cos(x))*exp(y+sin(x));});
            std::cout << "when n = " << n << ", 1-norm of error is: "
                << norm_1(test, [](double x, double y){return exp(y+sin(x));}) << std::endl;
            std::cout << "when n = " << n << ", 2-norm of error is: "
                << norm_2(test, [](double x, double y){return exp(y+sin(x));}) << std::endl;
            std::cout << "when n = " << n << ", inf-norm of error is: "
                << norm_inf(test, [](double x, double y){return exp(y+sin(x));}) << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << "! For The Irregular Domin:" << std::endl
              << std::endl;
    for(int i=0; i<3; i++){
        GET_THE_COEFFICIENT
        BoundaryCondi<domin::rmdisk> C(left_mixed(a1,b1),
                                       right_mixed(a2,b2),
                                       upper_mixed(a3,b3),
                                       lower_mixed(a4,b4),
                                       circle_mixed(a5,b5),
                                       j);
        disk D(2.0/3, 0.5, 0.25);
        std::vector<int> N = {8,16,32,64};
        std::cout << "No." << i+1 << std::endl;
        for(auto n: N){
            grid<domin::rmdisk> test(n,D);
            test.BVP_solver(C,[](double x, double y){return (sin(x)-1-cos(x)*cos(x))*exp(y+sin(x));});
            std::cout << "when n = " << n << ", 1-norm of error is: "
                << norm_1(test, [](double x, double y){return exp(y+sin(x));}) << std::endl;
            std::cout << "when n = " << n << ", 2-norm of error is: "
                << norm_2(test, [](double x, double y){return exp(y+sin(x));}) << std::endl;
            std::cout << "when n = " << n << ", inf-norm of error is: "
                << norm_inf(test, [](double x, double y){return exp(y+sin(x));}) << std::endl;
        }
        std::cout << std::endl;
    }
}

/**************************************************************
 > @description: The function is sin(x)+sin(y), 
 >               and the disk is (x-1)^2+(y-1/2)^2=(3/4)^2
**************************************************************/
void for_the_second_test()
{
    std::cout << "!! FOR THE 2ND TESTING FUNCTION: " << std::endl
              << "///////////////////////////////////////////////////"
              << std::endl << std::endl;
    
    auto left_mixed = [](double a, double b){
        return [a,b](double x, double y){
            return -cos(x)*a + b*(sin(x)+sin(y));
        };
    };
    auto right_mixed = [](double a, double b){
        return [a,b](double x, double y){
            return cos(x)*a + b*(sin(x)+sin(y));
        };
    };
    auto upper_mixed = [](double a, double b){
        return [a,b](double x, double y){
            return cos(y)*a + b*(sin(x)+sin(y));
        };
    };
    auto lower_mixed = [](double a, double b){
        return [a,b](double x, double y){
            return -cos(y)*a + b*(sin(x)+sin(y));
        };
    };
    auto circle_mixed = [](double a, double b){
        return [a,b](double x, double y){
            double x0 = 1, y0 = 0.5;
            double r = 0.75;
            return a*((x-x0)*cos(x)/r + (y-y0)*cos(y)/r) + b*(sin(x)+sin(y));
        };
    };

    std::cout << "! For The Regular Domin: " << std::endl 
              << std::endl;
    for(int i=0; i<3; i++){
        GET_THE_COEFFICIENT
        BoundaryCondi<domin::square> B(left_mixed(a1,b1),
                                       right_mixed(a2,b2),
                                       upper_mixed(a3,b3),
                                       lower_mixed(a4,b4),
                                       j);
        std::vector<int> N = {8,16,32,64};
        std::cout << "No." << i+1 << std::endl;
        for(auto n: N){
            grid<domin::square> test(n);
            test.BVP_solver(B, [](double x, double y){return sin(x)+sin(y);});
            std::cout << "when n = " << n << ", 1-norm of error is: "
                << norm_1(test, [](double x, double y){return sin(x)+sin(y);}) << std::endl;
            std::cout << "when n = " << n << ", 2-norm of error is: "
                << norm_2(test, [](double x, double y){return sin(x)+sin(y);}) << std::endl;
            std::cout << "when n = " << n << ", inf-norm of error is: "
                << norm_inf(test, [](double x, double y){return sin(x)+sin(y);}) << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << "! For The Irregular Domin:" << std::endl
              << std::endl;
    for(int i=0; i<3; i++){
        GET_THE_COEFFICIENT
        BoundaryCondi<domin::rmdisk> C(left_mixed(a1,b1),
                                       right_mixed(a2,b2),
                                       upper_mixed(a3,b3),
                                       lower_mixed(a4,b4),
                                       circle_mixed(a5,b5),
                                       j);
        disk D(1.0, 0.5, 0.75);
        std::vector<int> N = {8,16,32,64};
        std::cout << "No." << i+1 << std::endl;
        for(auto n: N){
            grid<domin::rmdisk> test(n,D);
            test.BVP_solver(C,[](double x, double y){return sin(x)+sin(y);});
            std::cout << "when n = " << n << ", 1-norm of error is: "
                << norm_1(test, [](double x, double y){return sin(x)+sin(y);}) << std::endl;
            std::cout << "when n = " << n << ", 2-norm of error is: "
                << norm_2(test, [](double x, double y){return sin(x)+sin(y);}) << std::endl;
            std::cout << "when n = " << n << ", inf-norm of error is: "
                << norm_inf(test, [](double x, double y){return sin(x)+sin(y);}) << std::endl;
        }
        std::cout << std::endl;
    }
}

/**************************************************************
 > @description: The function is x^2*y, 
 >               and the disk is (x-1/2)^2+(y-1/2)^2=(2/7)^2
**************************************************************/
void for_the_third_test()
{
    std::cout << "!! FOR THE 3RD TESTING FUNCTION: " << std::endl
              << "///////////////////////////////////////////////////"
              << std::endl << std::endl;
    
    auto left_mixed = [](double a, double b){
        return [a,b](double x, double y){
            return -2*x*y*a + x*x*y*b;
        };
    };
    auto right_mixed = [](double a, double b){
        return [a,b](double x, double y){
            return 2*x*y*a + x*x*y*b;
        };
    };
    auto upper_mixed = [](double a, double b){
        return [a,b](double x, double y){
            return x*x*a + x*x*y*b;
        };
    };
    auto lower_mixed = [](double a, double b){
        return [a,b](double x, double y){
            return -x*x*a + x*x*y*b;
        };
    };
    auto circle_mixed = [](double a, double b){
        return [a,b](double x, double y){
            double x0 = 0.5, y0 = 0.5;
            double r = 2.0/7;
            return a*((x-x0)*2*x*y/r + (y-y0)*x*x/r) + b*x*x*y;
        };
    };

    std::cout << "! For The Regular Domin: " << std::endl 
              << std::endl;
    for(int i=0; i<3; i++){
        GET_THE_COEFFICIENT
        BoundaryCondi<domin::square> B(left_mixed(a1,b1),
                                       right_mixed(a2,b2),
                                       upper_mixed(a3,b3),
                                       lower_mixed(a4,b4),
                                       j);
        std::vector<int> N = {8,16,32,64};
        std::cout << "No." << i+1 << std::endl;
        for(auto n: N){
            grid<domin::square> test(n);
            test.BVP_solver(B, [](double x, double y){return -2*y;});
            std::cout << "when n = " << n << ", 1-norm of error is: "
                << norm_1(test, [](double x, double y){return x*x*y;}) << std::endl;
            std::cout << "when n = " << n << ", 2-norm of error is: "
                << norm_2(test, [](double x, double y){return x*x*y;}) << std::endl;
            std::cout << "when n = " << n << ", inf-norm of error is: "
                << norm_inf(test, [](double x, double y){return x*x*y;}) << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << "! For The Irregular Domin:" << std::endl
              << std::endl;
    for(int i=0; i<3; i++){
        GET_THE_COEFFICIENT
        BoundaryCondi<domin::rmdisk> C(left_mixed(a1,b1),
                                       right_mixed(a2,b2),
                                       upper_mixed(a3,b3),
                                       lower_mixed(a4,b4),
                                       circle_mixed(a5,b5),
                                       j);
        disk D(0.5, 0.5, 2.0/7);
        std::vector<int> N = {8,16,32,64};
        std::cout << "No." << i+1 << std::endl;
        for(auto n: N){
            grid<domin::rmdisk> test(n,D);
            test.BVP_solver(C,[](double x, double y){return -2*y;});
            std::cout << "when n = " << n << ", 1-norm of error is: "
                << norm_1(test, [](double x, double y){return x*x*y;}) << std::endl;
            std::cout << "when n = " << n << ", 2-norm of error is: "
                << norm_2(test, [](double x, double y){return x*x*y;}) << std::endl;
            std::cout << "when n = " << n << ", inf-norm of error is: "
                << norm_inf(test, [](double x, double y){return x*x*y;}) << std::endl;
        }
        std::cout << std::endl;
    }
}

#undef GET_THE_COEFFICIENT

int main(int argc, char *argv[])
{
    for_the_first_test();
    for_the_second_test();
    for_the_third_test();
    return 0;
}