/**************************************************************
 > @Description    : one-dimensional multigrid method
 > @Version        : 1.0
 > @Author         : zhang-zh
 > @Date           : 2022-04-24 15:52
 > @LastEditTime   : 2022-04-26 21:31
**************************************************************/
#include<fstream>
#include"../include/GridDimOne.hpp"

using json = nlohmann::json;

/**************************************************************
 > @description: gradually reduce epsilon to 2.2e-16
**************************************************************/
void reduce_epsilon()
{
    std::cout << "////////////////////////////////////" << std::endl <<
        "!!! GRADUALLY REDUCE EPSILON !!!" << std::endl << std::endl;
    json j;
    std::ifstream jfile("./src/json/dim1BVC.json");
    jfile >> j;
    jfile.close();
    BoundaryCondition<1> B(-M_PI, 1, j);
    auto RightSideFunc = [](real x){
        return pow(M_PI,2)*exp(sin(M_PI*x))*(sin(M_PI*x)-pow(cos(M_PI*x),2));
    };
    auto U = [](real x){return exp(sin(M_PI*x));};

    Grid<1> grid1D(64,B,RightSideFunc);
    for(real epsilon=1e-8; epsilon>=1e-16; epsilon*=0.1){
        std::cout << "when epsilon = " << epsilon << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                            Operators::fullWeighting,Operators::linear,
                            100, epsilon);
        std::cout << std::endl;
        
    }

    std::cout << std::endl;
}

/**************************************************************
 > @description: f(x) = exp{sin(pi*x)}
**************************************************************/
void test_function_1()
{
    std::cout << "////////////////////////////////////" << std::endl <<
        "!!! FOR THE 1ST FUNCTION !!!" << std::endl << std::endl;
    json j;
    std::ifstream jfile("./src/json/dim1BVC.json");
    jfile >> j;
    jfile.close();
    BoundaryCondition<1> B(-M_PI, 1, j);
    auto RightSideFunc = [](real x){
        return pow(M_PI,2)*exp(sin(M_PI*x))*(sin(M_PI*x)-pow(cos(M_PI*x),2));
    };
    auto U = [](real x){return exp(sin(M_PI*x));};
#define DISPLAY_ERROR                                                                 \
    std::cout << "与解析解的误差的无穷范数为：" << grid1D.analytic_error(U) << std::endl;  \
    std::cout << "与离散精确解的误差的无穷范数为：" << grid1D.discrete_error() << std::endl;\
    std::cout << std::endl;
    for(int n=32; n<=256; n=n*2){
        std::cout << "!! When n = " << n << " !!" << std::endl;
        Grid<1> grid1D(n,B,RightSideFunc);
        std::cout << "! FMG, injection, linear !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::injection, Operators::linear);
        DISPLAY_ERROR
        std::cout << "! FMG, fullWeighting, linear !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::linear);
        DISPLAY_ERROR
        std::cout << "! FMG, injection, quadratic !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::injection, Operators::quadratic);
        DISPLAY_ERROR
        std::cout << "! FMG, fullWeighting, quadratic !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::quadratic);
        DISPLAY_ERROR
        std::cout << "! V-cycle, injection, linear !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::injection, Operators::linear,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, fullWeighting, linear !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::linear,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, injection, quadratic !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::injection, Operators::quadratic,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, fullWeighting, quadratic !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::quadratic,100,1e-8);
        DISPLAY_ERROR
        std::cout << std::endl;
    }
#undef DISPLAY_ERROR
    std::cout << "! error vector convergence rate !" << std::endl;
    Grid<1> grid(256,B,RightSideFunc);
    grid.error_convergence(U);
    std::cout << std::endl;
}

/**************************************************************
 > @description: f(x) = x*x-x
**************************************************************/
void test_function_2()
{
    std::cout << "////////////////////////////////////" << std::endl <<
        "!!! FOR THE 2ND FUNCTION !!!" << std::endl << std::endl;
    json j;
    std::ifstream jfile("./src/json/dim1BVC.json");
    jfile >> j;
    jfile.close();
    BoundaryCondition<1> B(1, 0, j);
    auto RightSideFunc = [](real x){
        return -2;
    };
    auto U = [](real x){return x*x-x;};
#define DISPLAY_ERROR                                                                 \
    std::cout << "与解析解的误差的无穷范数为：" << grid1D.analytic_error(U) << std::endl;  \
    std::cout << "与离散精确解的误差的无穷范数为：" << grid1D.discrete_error() << std::endl;\
    std::cout << std::endl;
    for(int n=32; n<=256; n=n*2){
        std::cout << "!! When n = " << n << " !!" << std::endl;
        Grid<1> grid1D(n,B,RightSideFunc);
        std::cout << "! FMG, injection, linear !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::injection, Operators::linear);
        DISPLAY_ERROR
        std::cout << "! FMG, fullWeighting, linear !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::linear);
        DISPLAY_ERROR
        std::cout << "! FMG, injection, quadratic !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::injection, Operators::quadratic);
        DISPLAY_ERROR
        std::cout << "! FMG, fullWeighting, quadratic !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::quadratic);
        DISPLAY_ERROR
        std::cout << "! V-cycle, injection, linear !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::injection, Operators::linear,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, fullWeighting, linear !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::linear,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, injection, quadratic !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::injection, Operators::quadratic,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, fullWeighting, quadratic !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::quadratic,100,1e-8);
        DISPLAY_ERROR
        std::cout << std::endl;
    }
#undef DISPLAY_ERROR
    std::cout << "! error vector convergence rate !" << std::endl;
    Grid<1> grid(256,B,RightSideFunc);
    grid.error_convergence(U);
    std::cout << std::endl;
}

/**************************************************************
 > @description: f(x) = sin(x)
**************************************************************/
void test_function_3()
{
    std::cout << "////////////////////////////////////" << std::endl <<
        "!!! FOR THE 3RD FUNCTION !!!" << std::endl << std::endl;
    json j;
    std::ifstream jfile("./src/json/dim1BVC.json");
    jfile >> j;
    jfile.close();
    BoundaryCondition<1> B(1, 0.5, j);
    auto RightSideFunc = [](real x){
        return -2.0/pow(x+1,3);
    };
    auto U = [](real x){return 1.0/(x+1);};
#define DISPLAY_ERROR                                                                 \
    std::cout << "与解析解的误差的无穷范数为：" << grid1D.analytic_error(U) << std::endl;  \
    std::cout << "与离散精确解的误差的无穷范数为：" << grid1D.discrete_error() << std::endl;\
    std::cout << std::endl;
    for(int n=32; n<=256; n=n*2){
        std::cout << "!! When n = " << n << " !!" << std::endl;
        Grid<1> grid1D(n,B,RightSideFunc);
        std::cout << "! FMG, injection, linear !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::injection, Operators::linear);
        DISPLAY_ERROR
        std::cout << "! FMG, fullWeighting, linear !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::linear);
        DISPLAY_ERROR
        std::cout << "! FMG, injection, quadratic !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::injection, Operators::quadratic);
        DISPLAY_ERROR
        std::cout << "! FMG, fullWeighting, quadratic !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::quadratic);
        DISPLAY_ERROR
        std::cout << "! V-cycle, injection, linear !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::injection, Operators::linear,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, fullWeighting, linear !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::linear,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, injection, quadratic !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::injection, Operators::quadratic,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, fullWeighting, quadratic !" << std::endl;
        grid1D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::quadratic,100,1e-8);
        DISPLAY_ERROR
        std::cout << std::endl;
    }
#undef DISPLAY_ERROR
    std::cout << "! error vector convergence rate !" << std::endl;
    Grid<1> grid(256,B,RightSideFunc);
    grid.error_convergence(U);
    std::cout << std::endl;
}

int main(int argc, char *argv[])
{
    reduce_epsilon();
    test_function_1();
    test_function_2();
    test_function_3();
    return 0;
}