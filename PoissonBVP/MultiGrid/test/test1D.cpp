#include<fstream>
#include"../include/GridDimOne.hpp"
#include"catch.hpp"

using json = nlohmann::json;

TEST_CASE("one-dimensional method", "[1D-multigrid]"){
    json j;
    std::ifstream jfile("./src/json/dim1BVC.json");
    jfile >> j;
    jfile.close();

    BoundaryCondition<1> B1(-M_PI, 1, j);
    auto RightSideFunc = [](real x){
        return pow(M_PI,2)*exp(sin(M_PI*x))*(sin(M_PI*x)-pow(cos(M_PI*x),2));
    };
    Grid<1> grid1D_1(64, B1, RightSideFunc);

    SECTION("test for V-cycle"){
        grid1D_1.Poisson_BVP_Multigrid(MultigridMethod::Vcycle,2.0/3,3,3,
                                   Operators::fullWeighting,Operators::linear,
                                   100, 1e-8);
        REQUIRE(grid1D_1.analytic_error([](real x){return exp(sin(M_PI*x));}) < 0.001);
    }

    SECTION("test for FMG"){
        grid1D_1.Poisson_BVP_Multigrid(MultigridMethod::FMG,2.0/3,3,3,
                                    Operators::fullWeighting, Operators::linear);
        REQUIRE(grid1D_1.analytic_error([](real x){return exp(sin(M_PI*x));}) < 0.1);
    }
}