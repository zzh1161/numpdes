/**************************************************************
 > @Description    : test for BVP
 > @Version        : 1.0
 > @Author         : zhang-zh
 > @Date           : 2022-03-30 16:29
 > @LastEditTime   : 2022-03-30 16:41
**************************************************************/
#include"../include/BVP.hpp"
#include"catch.hpp"
#include<fstream>

using namespace std;
using json = nlohmann::json;

TEST_CASE("BVP solver","[BVP_solver]"){

    json j;
    std::ifstream jfile("./test/testJSON/condition.json");
    jfile >> j;
    jfile.close();

    BoundaryCondi<domin::square> B([](double x, double y){return x*x+y*y;},
                                   [](double x, double y){return x*x+y*y;},
                                   [](double x, double y){return x*x+y*y;},
                                   [](double x, double y){return x*x+y*y;},
                                   j);

    disk D(1, 0.5, 0.75);
    BoundaryCondi<domin::rmdisk> C([](double x, double y){return x*x+y*y;},
                                   [](double x, double y){return x*x+y*y;},
                                   [](double x, double y){return x*x+y*y;},
                                   [](double x, double y){return x*x+y*y;},
                                   [](double x, double y){return x*x+y*y;},
                                   j);
    
    SECTION("test for regular domin"){
        grid<domin::square> G1(8);
        G1.BVP_solver(B, [](double x, double y){return -4;});
        REQUIRE(norm_inf(G1,[](double x, double y){return x*x+y*y;}) < 0.000000001);
    }

    SECTION("test for irregular domin"){
        grid<domin::rmdisk> G2(8,D);
        G2.BVP_solver(C, [](double x, double y){return -4;});
        REQUIRE(norm_inf(G2,[](double x, double y){return x*x+y*y;}) < 0.000000001);
    }

}