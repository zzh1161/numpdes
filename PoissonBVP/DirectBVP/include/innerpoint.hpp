/**************************************************************
 > @Description    : A header file for innerpoints
 > @Version        : 1.0
 > @Author         : zhang-zh
 > @Date           : 2022-03-06 20:24
 > @LastEditTime   : 2022-03-29 10:10
**************************************************************/
#ifndef __INNERPOINT_HPP__
#define __INNERPOINT_HPP__

#include<iostream>
#include<utility>
#include<cassert>
#include<cmath>

struct innerpoint
{
    enum positionType{
        corner     = 0, // points on the corner
        leftbound  = 1, // points on the left boundary
        rightbound = 2, // points on the right boundary
        upbound    = 3, // points on the upper boundary
        downbound  = 4, // points on the lower boundary
        ondisk     = 5, // points on the disk's boundary
        insquare   = 6, // points inside the square
        offsquare  = 7, // points off the square

        onXaxis    = 8, // points on a line paralleled to X-axis,
                        // with the line being part of the grid.
        onYaxis    = 9, // points on a line paralleled to Y-axis
                        // with the line being part of the grid.
        unknown    = 10
    };

            double _x,_y;  // coordinates of this point
    mutable double value;  // result after calculating, which means
                           // the unknown function's value at some inner point
    positionType   pos_type;

    innerpoint(){ value=0; }
    innerpoint(double x, double y): _x(x), _y(y){value=0;}
    innerpoint(double x, double y, positionType p): _x(x), _y(y), pos_type(p){value=0;}

    inline bool operator==(const innerpoint &t)const{
        return (abs(_x-t._x)+abs(_y-t._y) < 0.00000001);
    }

    double distance(const innerpoint &p)const{
        double dis_x = abs(_x-p._x);
        double dis_y = abs(_y-p._y);
        if(dis_x==0 || dis_y==0)
            return dis_x+dis_y;
        else
            return sqrt(pow(dis_x,2)+pow(dis_y,2));
    }
    
};


#endif