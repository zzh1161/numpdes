/**************************************************************
 > @Description    : output in a special way to plot
 > @Version        : 1.0
 > @Author         : zhang-zh
 > @Date           : 2022-04-25 20:17
 > @LastEditTime   : 2022-04-27 11:49
**************************************************************/
#include<fstream>
#include"../include/GridDimOne.hpp"
#include"../include/GridDimTwo.hpp"

using json = nlohmann::json;

/**************************************************************
 > @description: plot for one-dimensional multigrid
**************************************************************/
void one_dim_grid()
{
    std::cout << "/////////////// one-dimensional ///////////////" << std::endl;
    json j;
    std::ifstream jfile("./src/json/dim1BVC.json");
    jfile >> j;
    jfile.close();
    // 1st function: f(x) = exp{sin(pi*x)}
    std::cout << "1st function: " << std::endl;
    BoundaryCondition<1> B1(-M_PI, 1, j);
    auto RightSideFunc = [](real x){
        return pow(M_PI,2)*exp(sin(M_PI*x))*(sin(M_PI*x)-pow(cos(M_PI*x),2));
    };
    Grid<1> grid1D_1(256, B1, RightSideFunc);
    grid1D_1.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                   Operators::fullWeighting,Operators::linear,
                                   100, 1e-8);
    grid1D_1.display();
    std::cout << std::endl;
    // 2nd function: f(x) = x*x-x
    std::cout << "2nd function: " << std::endl;
    BoundaryCondition<1> B2(1, 0, j);
    Grid<1> grid1D_2(256, B2, [](real x){return -2;});
    grid1D_2.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                   Operators::fullWeighting,Operators::linear,
                                   100, 1e-8);
    grid1D_2.display();
    std::cout << std::endl;
    // 3rd function: f(x) = 1/(x+1)
    std::cout << "3rd function: " << std::endl;
    BoundaryCondition<1> B3(1, 0.5, j);
    Grid<1> grid1D_3(256, B3, [](real x){return -2.0/pow(x+1,3);});
    grid1D_3.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                   Operators::fullWeighting,Operators::linear,
                                   100, 1e-8);
    grid1D_3.display();
    std::cout << std::endl;
}

/**************************************************************
 > @description: plot for two-dimensional grid
**************************************************************/
void two_dim_grid()
{
    std::cout << "/////////////// two-dimensional ///////////////" << std::endl;
    json j;
    std::ifstream jfile("./src/json/dim2BVC.json");
    jfile >> j;
    jfile.close();
    // 1st function: f(x,y) = exp{y+sin(x)}
    std::cout << "1st function:" << std::endl;
    BoundaryCondition<2> B1([](real x, real y){return -cos(x)*exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           j);
    auto RightSideFunc1 = [](real x, real y){
        return exp(y+sin(x))*(sin(x)-1-pow(cos(x),2));
    };
    Grid<2> grid2D_1(64,B1,RightSideFunc1);
    grid2D_1.Poisson_BVP_Multigrid(MultigridMethod::Vcycle, 2.0/3,3,3,
                                Operators::fullWeighting, Operators::linear, 100,1e-8);
    grid2D_1.display();
    std::cout << std::endl;
    // 2nd function: f(x,y) = x^2+y^2
    std::cout << "2nd function: " << std::endl;
    BoundaryCondition<2> B2([](real x, real y){return -2*x;},
                           [](real x, real y){return x*x+y*y;},
                           [](real x, real y){return x*x+y*y;},
                           [](real x, real y){return x*x+y*y;},
                           j);
    auto RightSideFunc2 = [](real x, real y){
        return -4;
    };
    Grid<2> grid2D_2(64,B2,RightSideFunc2);
    grid2D_2.Poisson_BVP_Multigrid(MultigridMethod::Vcycle, 2.0/3,3,3,
                                Operators::fullWeighting, Operators::linear, 100,1e-8);
    grid2D_2.display();
    std::cout << std::endl;
    // 3rd function: f(x,y) = sin(x+y)
    std::cout << "3rd function:" << std::endl;
    BoundaryCondition<2> B3([](real x, real y){return -cos(x+y);},
                           [](real x, real y){return sin(x+y);},
                           [](real x, real y){return sin(x+y);},
                           [](real x, real y){return sin(x+y);},
                           j);
    auto RightSideFunc3 = [](real x, real y){
        return 2*sin(x+y);
    };
    Grid<2> grid2D_3(64,B3,RightSideFunc3);
    grid2D_3.Poisson_BVP_Multigrid(MultigridMethod::Vcycle, 2.0/3,3,3,
                                Operators::fullWeighting, Operators::linear, 100,1e-8);
    grid2D_3.display();
    std::cout << std::endl;
}

int main(int argc, char *argv[])
{
    one_dim_grid();
    two_dim_grid();
    return 0;
}