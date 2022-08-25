/**************************************************************
 > @Description    : test for Boundary conditions
 > @Version        : 1.0
 > @Author         : zhang-zh
 > @Date           : 2022-03-30 16:16
 > @LastEditTime   : 2022-03-30 16:28
**************************************************************/
#include"../include/boundcondi.hpp"
#include"catch.hpp"
#include<fstream>

using namespace std;
using json = nlohmann::json;

TEST_CASE("Generated BoundaryCondi with JSON","[Boundary_condition]"){

    json j;
    std::ifstream jfile("./test/testJSON/condition.json");
    jfile >> j;
    jfile.close();

    BoundaryCondi<domin::square> B1([](double x, double y){return x+y;},
                                    [](double x, double y){return x-y;},
                                    [](double x, double y){return x*y;},
                                    [](double x, double y){return x/y;},
                                    j);
    BoundaryCondi<domin::rmdisk> B2([](double x, double y){return x+y;},
                                    [](double x, double y){return x-y;},
                                    [](double x, double y){return x*y;},
                                    [](double x, double y){return x/y;},
                                    [](double x, double y){return x*y;},
                                    j);
    SECTION("whether function is correct"){
        REQUIRE(B1.Lboundary.first(1,2) == 3);
        REQUIRE(B2.Lboundary.first(2,3) == 5);
    }
}