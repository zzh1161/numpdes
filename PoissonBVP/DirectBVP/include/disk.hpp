/**************************************************************
 > @Description    : A header file for disk
 > @Version        : 1.0
 > @Author         : zhang-zh
 > @Date           : 2022-03-06 21:32
 > @LastEditTime   : 2022-03-27 09:47
**************************************************************/
#ifndef __DISK_HPP__
#define __DISK_HPP__

#include<iostream>
#include<cassert>
#include<optional> // C++17
#include<cmath>

class disk
{
public:
    disk():_x0(0), _y0(0), _r(0){}
    disk(double x, double y, double r):_x0(x), _y0(y), _r(r){}

    inline bool whether_in(double x, double y) const {
        return (x-_x0)*(x-_x0)+(y-_y0)*(y-_y0) <= _r*_r;
    }
    inline bool whether_out(double x, double y) const {
        return (x-_x0)*(x-_x0)+(y-_y0)*(y-_y0) > _r*_r;
    }
    inline bool whether_onit(double x, double y) const {
        return abs((x-_x0)*(x-_x0)+(y-_y0)*(y-_y0)-_r*_r)<0.00000001;
    }

    std::optional<std::pair<double,double>> if_cross_axis(
                double x1, double y1, double x2, double y2) const {
        if(x1 != x2) assert(y1 == y2);
        else         assert(y1 != y2);
        if(x1 != x2){
            if(_r*_r <= (y1-_y0)*(y1-_y0))
                return std::nullopt;
            else{
                double e1 = _x0 + sqrt(pow(_r,2)-pow(y1-_y0,2));
                double e2 = _x0 - sqrt(pow(_r,2)-pow(y1-_y0,2));
                if((x1-e1)*(x2-e1)>0 && (x1-e2)*(x2-e2)>0)
                    return std::nullopt;
                else if((x1-e1)*(x2-e1)<=0 && (x1-e2)*(x2-e2)>0)
                    return std::pair<double,double>{e1,y1};
                else if((x1-e1)*(x2-e1)>0 && (x1-e2)*(x2-e2)<=0)
                    return std::pair<double,double>{e2,y1};
                else
                    return abs(x1-e1)<abs(x1-e2)?
                           std::pair<double,double>{e1,y1}:
                           std::pair<double,double>{e2,y1};
            }
        }else{
            if(pow(_r,2) <= pow(x1-_x0,2))
                return std::nullopt;
            else{
                double e1 = _y0 + sqrt(pow(_r,2)-pow(x1-_x0,2));
                double e2 = _y0 - sqrt(pow(_r,2)-pow(x1-_x0,2));
                if((y1-e1)*(y2-e1)>0 && (y1-e2)*(y2-e2)>0)
                    return std::nullopt;
                else if((y1-e1)*(y2-e1)<=0 && (y1-e2)*(y2-e2)>0)
                    return std::pair<double,double>{x1,e1};
                else if((y1-e1)*(y2-e1)>0 && (y1-e2)*(y2-e2)<=0)
                    return std::pair<double,double>{x1,e2};
                else 
                    return abs(y1-e1)<abs(y1-e2)?
                           std::pair<double,double>{x1,e1}:
                           std::pair<double,double>{x1,e2};
            }
        }
    }

// private:
    double _x0,_y0;  // center coordinates
    double _r;       // radius
};


#endif