#include<fstream>
#include"../include/GridDimTwo.hpp"
#include"catch.hpp"

using json = nlohmann::json;

TEST_CASE("two-dimensional method", "[2D-multigrid]"){
    json j;
    std::ifstream jfile("./src/json/dim2BVC.json");
    jfile >> j;
    jfile.close();

    BoundaryCondition<2> B1([](real x, real y){return -cos(x)*exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           [](real x, real y){return exp(y+sin(x));},
                           j);
    auto RightSideFunc1 = [](real x, real y){
        return exp(y+sin(x))*(sin(x)-1-pow(cos(x),2));
    };
    Grid<2> grid2D_1(64,B1,RightSideFunc1);

    SECTION("test for V-cycle"){
        grid2D_1.Poisson_BVP_Multigrid(MultigridMethod::Vcycle, 2.0/3,3,3,
                                Operators::fullWeighting, Operators::linear, 100,1e-8);
        REQUIRE(grid2D_1.analytic_error([](real x, real y){return exp(y+sin(x));}) < 0.0001);
    }

    SECTION("test for FMG"){
        grid2D_1.Poisson_BVP_Multigrid(MultigridMethod::FMG, 2.0/3,3,3,
                                Operators::fullWeighting, Operators::linear);
        REQUIRE(grid2D_1.analytic_error([](real x, real y){return exp(y+sin(x));}) < 0.1);
    }
}