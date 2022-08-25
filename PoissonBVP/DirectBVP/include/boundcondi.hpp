/**************************************************************
 > @Description    : boundary conditions
 > @Version        : 1.0
 > @Author         : zhang-zh
 > @Date           : 2022-03-13 11:05
 > @LastEditTime   : 2022-03-29 22:11
**************************************************************/
#ifndef __BOUNDARY_CONDITIONS_HPP__
#define __BOUNDARY_CONDITIONS_HPP__

#include<iostream>
#include<string>
#include<utility>
#include<vector>
#include<functional>
#include"domin.hpp"
#include"json.hpp"

enum class boundtype{
    Dirichlet = 0,
    Neumann   = 1 
};

template<domin condi>
class BoundaryCondi{};

template<>
class BoundaryCondi<domin::square>
{
public:
    using json = nlohmann::json;

    BoundaryCondi(){}
    BoundaryCondi(auto L, auto R, auto U, auto D, json &j){
        Lboundary.first = L; Lboundary.second = j.at("BoundaryOfSquare").at("left").at("BoundaryType");
        Rboundary.first = R; Rboundary.second = j.at("BoundaryOfSquare").at("right").at("BoundaryType");
        Uboundary.first = U; Uboundary.second = j.at("BoundaryOfSquare").at("up").at("BoundaryType");
        Dboundary.first = D; Dboundary.second = j.at("BoundaryOfSquare").at("down").at("BoundaryType");

        LMix.first  = j.at("BoundaryOfSquare").at("left").at("CoeffiOfNeum");
        LMix.second = j.at("BoundaryOfSquare").at("left").at("CoeffiOfDirch");

        RMix.first  = j.at("BoundaryOfSquare").at("right").at("CoeffiOfNeum");
        RMix.second = j.at("BoundaryOfSquare").at("right").at("CoeffiOfDirch");

        UMix.first  = j.at("BoundaryOfSquare").at("up").at("CoeffiOfNeum");
        UMix.second = j.at("BoundaryOfSquare").at("up").at("CoeffiOfDirch");

        DMix.first  = j.at("BoundaryOfSquare").at("down").at("CoeffiOfNeum");
        DMix.second = j.at("BoundaryOfSquare").at("down").at("CoeffiOfDirch");

        value_at_00 = j.at("ValueAt00");
    }

    std::pair<std::function<double(double,double)>,boundtype> 
                            Lboundary, Rboundary, Uboundary, Dboundary;
    // coefiicients for the third situation: Mixed
    std::pair<double,double> LMix, RMix, UMix, DMix;
    // if they all are Neumann's, we have to give a known point.
    double value_at_00;
};

template<>
class BoundaryCondi<domin::rmdisk>
{
public:
    using json = nlohmann::json;

    BoundaryCondi(){}
    BoundaryCondi(auto L, auto R, auto U, auto D, auto C, json &j){
        Lboundary.first = L; Lboundary.second = j.at("BoundaryOfSquare").at("left").at("BoundaryType");
        Rboundary.first = R; Rboundary.second = j.at("BoundaryOfSquare").at("right").at("BoundaryType");
        Uboundary.first = U; Uboundary.second = j.at("BoundaryOfSquare").at("up").at("BoundaryType");
        Dboundary.first = D; Dboundary.second = j.at("BoundaryOfSquare").at("down").at("BoundaryType");
        Cirboundary.first = C; Cirboundary.second = j.at("BoundaryOfDisk").at("BoundaryType");

        LMix.first  = j.at("BoundaryOfSquare").at("left").at("CoeffiOfNeum");
        LMix.second = j.at("BoundaryOfSquare").at("left").at("CoeffiOfDirch");

        RMix.first  = j.at("BoundaryOfSquare").at("right").at("CoeffiOfNeum");
        RMix.second = j.at("BoundaryOfSquare").at("right").at("CoeffiOfDirch");

        UMix.first  = j.at("BoundaryOfSquare").at("up").at("CoeffiOfNeum");
        UMix.second = j.at("BoundaryOfSquare").at("up").at("CoeffiOfDirch");

        DMix.first  = j.at("BoundaryOfSquare").at("down").at("CoeffiOfNeum");
        DMix.second = j.at("BoundaryOfSquare").at("down").at("CoeffiOfDirch");

        CMix.first  = j.at("BoundaryOfDisk").at("CoeffiOfNeum");
        CMix.second = j.at("BoundaryOfDisk").at("CoeffiOfDirch");
    
        value_at_00 = j.at("ValueAt00");
    }

    std::pair<std::function<double(double,double)>,boundtype> 
                        Lboundary, Rboundary, Uboundary, Dboundary, Cirboundary;
    // coefiicients for the third situation: Mixed
    std::pair<double,double> LMix, RMix, UMix, DMix, CMix;
    // if they all are Neumann's, we have to give a known point.
    double value_at_00;
};


#endif