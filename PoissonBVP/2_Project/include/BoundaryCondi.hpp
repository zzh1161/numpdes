/**************************************************************
 > @Description    : Boundary Conditions
 > @Version        : 1.0
 > @Author         : zhang-zh
 > @Date           : 2022-04-12 21:57
 > @LastEditTime   : 2022-04-23 17:19
**************************************************************/
#ifndef __PRO2_BOUNDARYCONDI_HPP__
#define __PRO2_BOUNDARYCONDI_HPP__

#include<algorithm>
#include<functional>
#include<iostream>
#include<memory>
#include<numeric>
#include<utility>
#include"BeforeAll.hpp"
#include"json.hpp"

template<int dim, Domin domin = Domin::regular>
struct BoundaryCondition{};

template<>
struct BoundaryCondition<1>
{
public:
    BoundaryCondition(){}
    BoundaryCondition(real l, real r, nlohmann::json &j){
        leftMixedCoeffi.first   = j["BoundForOneDim"]["Left"]["CoeffiOfNeumann"];
        leftMixedCoeffi.second  = j["BoundForOneDim"]["Left"]["CoeffiOfDirichlet"];
        rightMixedCoeffi.first  = j["BoundForOneDim"]["Right"]["CoeffiOfNeumann"];
        rightMixedCoeffi.second = j["BoundForOneDim"]["Right"]["CoeffiOfDirichlet"];
        leftBound.first = l; rightBound.first = r;
        leftBound.second  = j["BoundForOneDim"]["Left"]["Type"];
        rightBound.second = j["BoundForOneDim"]["Right"]["Type"];
    }

    std::pair<real,real> leftMixedCoeffi, rightMixedCoeffi;
    std::pair<real,BoundCondiType> leftBound, rightBound;
};

template<>
struct BoundaryCondition<2>
{
    BoundaryCondition(){}
    BoundaryCondition(auto l, auto r, auto u, auto d, nlohmann::json &j){
        leftMixed.first   = j["BoundForTwoDim"]["Left"]["CoeffiOfNeumann"];
        leftMixed.second  = j["BoundForTwoDim"]["Left"]["CoeffiOfDirichlet"];
        rightMixed.first  = j["BoundForTwoDim"]["Right"]["CoeffiOfNeumann"];
        rightMixed.second = j["BoundForTwoDim"]["Right"]["CoeffiOfDirichlet"];
        upperMixed.first  = j["BoundForTwoDim"]["Upper"]["CoeffiOfNeumann"];
        upperMixed.second = j["BoundForTwoDim"]["Upper"]["CoeffiOfDirichlet"];
        downMixed.first   = j["BoundForTwoDim"]["Down"]["CoeffiOfNeumann"];
        downMixed.second  = j["BoundForTwoDim"]["Down"]["CoeffiOfDirichlet"];
        leftBound  = l;
        rightBound = r;
        upperBound = u;
        downBound  = d;
    }
    
    std::function<real(real,real)> leftBound, rightBound, upperBound, downBound;
    std::pair<real,real>           leftMixed, rightMixed, upperMixed, downMixed;
};



#endif