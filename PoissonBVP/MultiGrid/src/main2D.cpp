/**************************************************************
 > @Description    : two dimensional multigrid method
 > @Version        : 1.0
 > @Author         : zhang-zh
 > @Date           : 2022-04-24 22:05
 > @LastEditTime   : 2022-04-27 15:12
**************************************************************/
#include<chrono>
#include<ctime>
#include<fstream>
#include"../include/GridDimTwo.hpp"

using json = nlohmann::json;

/**************************************************************
 > @description: compare running time between two methods
**************************************************************/
void CPU_time()
{
    clock_t begin, end;
    std::cout << "////////////////////////////////////" << std::endl <<
        "!!! WHEN n=128, CPU TIME COMPARISON !!!" << std::endl << std::endl;
    json j;
    std::ifstream jfile("./src/json/dim2BVC.json");
    jfile >> j;
    jfile.close();
    BoundaryCondition<2> B([](real x, real y){return -cos(x)*exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           j);
    auto RightSideFunc = [](real x, real y){
        return exp(y+sin(x))*(sin(x)-1-pow(cos(x),2));
    };
    auto U = [](real x, real y){
        return exp(y+sin(x));
    };

    Grid<2> grid2D(128, B, RightSideFunc);

    begin = clock();
    grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                Operators::injection, Operators::linear,100,1e-8);
    end = clock();
    std::cout << "多重网格法用时：" << (end-begin)*1.0/CLOCKS_PER_SEC << std::endl;

    begin = clock();
    std::cout << "与离散准确解的误差: " << grid2D.discrete_error() << std::endl;
    end = clock();
    std::cout << "直接求解法用时：" << (end-begin)*1.0/CLOCKS_PER_SEC << std::endl;

    std::cout << std::endl;
}

/**************************************************************
 > @description: gradually reduce epsilon to 2.2e-16
**************************************************************/
void reduce_epsilon()
{
    std::cout << "////////////////////////////////////" << std::endl <<
        "!!! GRADUALLY REDUCE EPSILON !!!" << std::endl << std::endl;
    json j;
    std::ifstream jfile("./src/json/dim2BVC.json");
    jfile >> j;
    jfile.close();
    BoundaryCondition<2> B([](real x, real y){return -cos(x)*exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           j);
    auto RightSideFunc = [](real x, real y){
        return exp(y+sin(x))*(sin(x)-1-pow(cos(x),2));
    };
    auto U = [](real x, real y){
        return exp(y+sin(x));
    };

    Grid<2> grid2D(64,B,RightSideFunc);
    for(real epsilon=1e-8; epsilon>=1e-16; epsilon*=0.1){
        std::cout << "when epsilon = " << epsilon << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                            Operators::fullWeighting,Operators::linear,
                            100, epsilon);
        std::cout << "与解析解的误差的无穷范数为: " << grid2D.analytic_error(U) << std::endl;
        std::cout << std::endl;
        
    }
    std::cout << std::endl;
}

/**************************************************************
 > @description: f(x,y) = exp{y+sin(x)}
**************************************************************/
void test_function_1()
{
    std::cout << "////////////////////////////////////" << std::endl <<
        "!!! FOR THE 1ST FUNCTION !!!" << std::endl << std::endl;
    json j;
    std::ifstream jfile("./src/json/dim2BVC.json");
    jfile >> j;
    jfile.close();
    BoundaryCondition<2> B([](real x, real y){return -cos(x)*exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           j);
    auto RightSideFunc = [](real x, real y){
        return exp(y+sin(x))*(sin(x)-1-pow(cos(x),2));
    };
    auto U = [](real x, real y){
        return exp(y+sin(x));
    };
#define DISPLAY_ERROR                                                                \
    std::cout << "与解析解的误差的无穷范数为：" << grid2D.analytic_error(U) << std::endl; \
    std::cout << std::endl;
    // std::cout << "与离散精确解的相对误差为：" << grid2D.discrete_error() << std::endl;
    for(int n=32; n<=256; n*=2){
        std::cout << "!! When n = " << n << " !!" << std::endl;
        Grid<2> grid2D(n, B, RightSideFunc);
        std::cout << "! FMG, injection, linear !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::injection, Operators::linear);
        DISPLAY_ERROR
        std::cout << "! FMG, fullWeighting, linear !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::linear);
        DISPLAY_ERROR
        std::cout << "! FMG, injection, quadratic !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::injection, Operators::quadratic);
        DISPLAY_ERROR
        std::cout << "! FMG, fullWeighting, quadratic !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::quadratic);
        DISPLAY_ERROR
        std::cout << "! V-cycle, injection, linear !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::injection, Operators::linear,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, fullWeighting, linear !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::linear,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, injection, quadratic !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::injection, Operators::quadratic,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, fullWeighting, quadratic !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::quadratic,100,1e-8);
        DISPLAY_ERROR
        std::cout << std::endl;
    }
#undef DISPLAY_ERROR
    std::cout << "! error vector convergence rate !" << std::endl;
    Grid<2> grid(256,B,RightSideFunc);
    grid.error_convergece(U);
    std::cout << std::endl;
}

/**************************************************************
 > @description: f(x,y) = x^2+y^2
**************************************************************/
void test_function_2()
{
    std::cout << "////////////////////////////////////" << std::endl <<
        "!!! FOR THE 2ND FUNTION !!!" << std::endl << std::endl;
    json j;
    std::ifstream jfile("./src/json/dim2BVC.json");
    jfile >> j;
    jfile.close();
    BoundaryCondition<2> B([](real x, real y){return -2*x;},
                           [](real x, real y){return x*x+y*y;},
                           [](real x, real y){return x*x+y*y;},
                           [](real x, real y){return x*x+y*y;},
                           j);
    auto RightSideFunc = [](real x, real y){
        return -4;
    };
    auto U = [](real x, real y){
        return x*x+y*y;
    };
#define DISPLAY_ERROR                                                                \
    std::cout << "与解析解的误差的无穷范数为：" << grid2D.analytic_error(U) << std::endl; \
    std::cout << std::endl; 
    // std::cout << "与离散精确解的相对误差为：" << grid2D.discrete_error() << std::endl;
    for(int n=32; n<=256; n*=2){
        std::cout << "!! When n = " << n << " !!" << std::endl;
        Grid<2> grid2D(n, B, RightSideFunc);
        std::cout << "! FMG, injection, linear !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::injection, Operators::linear);
        DISPLAY_ERROR
        std::cout << "! FMG, fullWeighting, linear !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::linear);
        DISPLAY_ERROR
        std::cout << "! FMG, injection, quadratic !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::injection, Operators::quadratic);
        DISPLAY_ERROR
        std::cout << "! FMG, fullWeighting, quadratic !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::quadratic);
        DISPLAY_ERROR
        std::cout << "! V-cycle, injection, linear !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::injection, Operators::linear,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, fullWeighting, linear !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::linear,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, injection, quadratic !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::injection, Operators::quadratic,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, fullWeighting, quadratic !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::quadratic,100,1e-8);
        DISPLAY_ERROR
        std::cout << std::endl;
    }
#undef DISPLAY_ERROR
    std::cout << "! error vector convergence rate !" << std::endl;
    Grid<2> grid(256,B,RightSideFunc);
    grid.error_convergece(U);
    std::cout << std::endl;
}

/**************************************************************
 > @description: f(x,y) = sin(x+y)
**************************************************************/
void test_function_3()
{
    std::cout << "////////////////////////////////////" << std::endl <<
        "!!! FOR THE 3RD FUNTION !!!" << std::endl << std::endl;
    json j;
    std::ifstream jfile("./src/json/dim2BVC.json");
    jfile >> j;
    jfile.close();
    BoundaryCondition<2> B([](real x, real y){return -cos(x+y);},
                           [](real x, real y){return sin(x+y);},
                           [](real x, real y){return sin(x+y);},
                           [](real x, real y){return sin(x+y);},
                           j);
    auto RightSideFunc = [](real x, real y){
        return 2*sin(x+y);
    };
    auto U = [](real x, real y){
        return sin(x+y);
    };
#define DISPLAY_ERROR                                                                \
    std::cout << "与解析解的误差的无穷范数为：" << grid2D.analytic_error(U) << std::endl; \
    std::cout << std::endl;
    // std::cout << "与离散精确解的相对误差为：" << grid2D.discrete_error() << std::endl;
    for(int n=32; n<=256; n*=2){
        std::cout << "!! When n = " << n << " !!" << std::endl;
        Grid<2> grid2D(n, B, RightSideFunc);
        std::cout << "! FMG, injection, linear !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::injection, Operators::linear);
        DISPLAY_ERROR
        std::cout << "! FMG, fullWeighting, linear !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::linear);
        DISPLAY_ERROR
        std::cout << "! FMG, injection, quadratic !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::injection, Operators::quadratic);
        DISPLAY_ERROR
        std::cout << "! FMG, fullWeighting, quadratic !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::quadratic);
        DISPLAY_ERROR
        std::cout << "! V-cycle, injection, linear !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::injection, Operators::linear,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, fullWeighting, linear !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::linear,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, injection, quadratic !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::injection, Operators::quadratic,100,1e-8);
        DISPLAY_ERROR
        std::cout << "! V-cycle, fullWeighting, quadratic !" << std::endl;
        grid2D.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::quadratic,100,1e-8);
        DISPLAY_ERROR
        std::cout << std::endl;
    }
#undef DISPLAY_ERROR
    std::cout << "! error vector convergence rate !" << std::endl;
    Grid<2> grid(256,B,RightSideFunc);
    grid.error_convergece(U);
    std::cout << std::endl;
}

int main(int argc, char *agrv[])
{
    CPU_time();
    reduce_epsilon();
    test_function_1();
    test_function_2();
    test_function_3();
    return 0;
}