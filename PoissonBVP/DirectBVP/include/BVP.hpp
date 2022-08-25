/**************************************************************
 > @Description    : class grid and BVP methods
 > @Version        : 1.0
 > @Author         : zhang-zh
 > @Date           : 2022-03-06 20:34
 > @LastEditTime   : 2022-03-30 14:36
**************************************************************/
#ifndef __GRID_BVP_HPP__
#define __GRID_BVP_HPP__

#include<iostream>
#include<vector>
#include<utility>
#include<cassert>
#include<algorithm>
#include<numeric>
#include<map>
#include<lapacke.h>
#include<initializer_list>
#include"innerpoint.hpp"
#include"disk.hpp"
#include"boundcondi.hpp"

/**
 * @tparam size   : size n (n*n grid)
 * @tparam condi  : which kind of domins
 */
template<domin condi>
class grid{};

//////////////////////////////////////////////////////////////////////////////

template<>
class grid<domin::square>
{
public:
    grid(int _size){
        size = _size;
        assert(size > 2);
        int n = size;
        double h = 1.0/(n-1);
        points.resize(n*n);
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                points[i*n+j]._x = j*h;
                points[i*n+j]._y = i*h;
            }
        }
    }
    
    /**************************************************************
     > @description: BVP solver for regular domin
    **************************************************************/    
    void BVP_solver(BoundaryCondi<domin::square> BdC,
                    std::function<double(double,double)> f);
    
    friend std::ostream& operator<<(std::ostream& os,
                            grid<domin::square> &g){
        int n = g.size;
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                os << g.points[(n-1-i)*n+j].value << ", ";
            }
            os << std::endl;
        }
        return os;
    }

    friend double norm_1(const auto &g, auto func);
    friend double norm_2(const auto &g, auto func);
    friend double norm_inf(const auto &g, auto func);
    
private:
    int size;
    mutable std::vector<innerpoint> points;
};

//////////////////////////////////////////////////////////////////////////////

template<>
class grid<domin::rmdisk>
{
public:
    using coordinates = std::pair<double,double>;

    grid(int _size, disk d)
    {
        size = _size;
        assert(size > 2);
        D = d;
        int     n = size;
        double  h = 1.0/(n-1);
        int count = 0; // a flag
        
        for(int j=0; j<n; j++){
            for(int i=0; i<n; i++){
                double x=i*h, y=j*h;
                if(D.whether_out(x,y)){
                    innerpoint neigh[5] = {
                        innerpoint(x,y,type_of_point(x,y)),
                        left_pnt(x,y),
                        right_pnt(x,y),
                        upper_pnt(x,y),
                        down_pnt(x,y)
                    };
                    for(int k=0; k<5; k++){
                        if(neigh[k].pos_type != innerpoint::offsquare){
                            for(auto item: points){
                                if(item == neigh[k]){
                                    count=1; break;
                                }
                            }
                            if(count != 1){
                                points.push_back(neigh[k]);
                            }
                            count = 0;
                        }
                    }
                }
            }
        }
    }

    void display()const{
        // const char *postypechar[] = {
        //     "corner", "leftbound", "rightbound", "upbound",
        //     "downbound", "ondisk", "insquare", "offsquare"
        // };
        // for(int i=0; i<points.size(); i++){
        //     std::cout << i << ": (" << points[i]._x << ", " <<
        //         points[i]._y << ") " << postypechar[points[i].pos_type] 
        //         << std::endl; 
        // }
        std::cout << "x = [";
        for(auto p : points)
            std::cout << p._x << ", ";
        std::cout << "]" << std::endl << "y = [";
        for(auto p : points)
            std::cout << p._y << ", ";
        std::cout << "]" << std::endl << "z = [";
        for(auto p : points)
            std::cout << p.value << ", ";
        std::cout << "]" << std::endl;
    }

    friend std::ostream& operator<<(std::ostream& os,
                                const grid<domin::rmdisk> &g){
        for(auto point : g.points){
            os << "(" << point._x << ", " << point._y
               << "): " << point.value << std::endl;
        }
        return os;
    }

    /**************************************************************
     > @description: BVP solver for irregular domin
    **************************************************************/    
    void BVP_solver(BoundaryCondi<domin::rmdisk> BdC,
                    std::function<double(double,double)> f);

    friend double norm_1(const auto &g, auto func);
    friend double norm_2(const auto &g, auto func);
    friend double norm_inf(const auto &g, auto func);

private:
    int size;
    mutable std::vector<innerpoint> points;
    disk D;

    int ret_order(innerpoint &p)const{
        for(int i=0; i<points.size(); i++){
            if(p == points[i])
                return i;
        }
        return -1;
    }

    innerpoint find(double x, double y)const{
        innerpoint temp(x,y,innerpoint::unknown);
        for(auto i : points)
            if(i == temp)
                return i;
        return temp;
    }

    innerpoint intersection_with_grid(double x1, double y1,
                                      double x2, double y2){
        // with (x1,y1) center of the disk,
        // and (x2,y2) a point on the disk.
        double h = 1.0/(size-1);
        int n = size;
        std::vector<std::pair<double,innerpoint::positionType>> many_k;
        for(int i=0; i<n; i++){
            double k = (i*h-x1)/(x2-x1);
            if(k > 1)
                many_k.push_back({k,innerpoint::onYaxis});
        }
        for(int i=0; i<n; i++){
            double k = (i*h-y1)/(y2-y1);
            if(k > 1)
                many_k.push_back({k,innerpoint::onXaxis});
        }
        auto min_k = many_k[0];
        for(auto i : many_k){
            if(i.first < min_k.first)
                min_k = i;
        }
        double x3 = min_k.first*(x2-x1)+x1;
        double y3 = min_k.first*(y2-y1)+y1;
        return innerpoint(x3,y3,min_k.second);
    }

    innerpoint::positionType type_of_point(double x, double y)const{
        if(D.whether_onit(x,y))
            return innerpoint::ondisk;
        else{
            if((x==0&&y==0)||(x==0&&y==1)||(x==1&&y==0)||(x==1&&y==1))
                return innerpoint::corner;
            else if(x==0)
                return innerpoint::leftbound;
            else if(x==1)
                return innerpoint::rightbound;
            else if(y==0)
                return innerpoint::downbound;
            else if(y==1)
                return innerpoint::upbound;
            else if(x<0 || x>1 || y<0 || y>1)
                return innerpoint::offsquare;
            else
                return innerpoint::insquare;
        }
    }

#define OPERATION_RET(a,b)                          \
    assert(D.whether_out(x,y));                     \
    auto res = D.if_cross_axis(x,y,a,b);            \
    if(res.has_value())                             \
        return innerpoint(res->first,res->second,   \
                               innerpoint::ondisk); \
    else                                            \
        return innerpoint(a,b,type_of_point(a,b));
    // assume arguments are the points on the grid
    innerpoint left_pnt(double x, double y)const{
        double h = 1.0/(size-1);
        if(x >= h){
            OPERATION_RET(x-h,y)
        }else
            return innerpoint(x-h,y,innerpoint::offsquare);
    }
    innerpoint right_pnt(double x, double y)const{
        double h = 1.0/(size-1);
        if(x+h <= 1){
            OPERATION_RET(x+h,y)
        }else
            return innerpoint(x+h,y,innerpoint::offsquare);
    }
    innerpoint upper_pnt(double x, double y)const{
        double h = 1.0/(size-1);
        if(y+h <= 1){
            OPERATION_RET(x,y+h)
        }else
            return innerpoint(x,y+h,innerpoint::offsquare);
    }
    innerpoint down_pnt(double x, double y)const{
        double h = 1.0/(size-1);
        if(y >= h){
            OPERATION_RET(x,y-h)
        }else
            return innerpoint(x,y-h,innerpoint::offsquare);
    }
#undef OPERATION_RET
};

////////////////////////////////////////////////////////////////////////////////
/***************************** norm function **********************************/
////////////////////////////////////////////////////////////////////////////////

#define GENERATE_ERROR                         \
    int n = g.size;                            \
    double h = 1.0/(n-1);                      \
    std::vector<double> err;                   \
    for(auto i : g.points){                    \
        double x = i._x;                       \
        double y = i._y;                       \
        err.push_back(abs(i.value-func(x,y))); \
    }

double norm_1(const auto &g, auto func)
{
    GENERATE_ERROR
    return h*std::accumulate(
                    err.begin(),
                    err.end(),
                    0.0
                    );
}

double norm_2(const auto &g, auto func)
{
    GENERATE_ERROR
    std::transform(
            err.begin(),
            err.end(),
            err.begin(),
            [](double x){return x*x;}
            );
    return sqrt(h*std::accumulate(
                    err.begin(),
                    err.end(),
                    0.0
                    ));
}

double norm_inf(const auto &g, auto func)
{
    GENERATE_ERROR
    auto MAX = std::max_element(err.begin(), err.end());
    return *MAX;
}

#undef GENERATE_ERROR

////////////////////////////////////////////////////////////////////////////////
/******************************* BVP solver ***********************************/
////////////////////////////////////////////////////////////////////////////////

void grid<domin::square>::BVP_solver(
                    BoundaryCondi<domin::square> BdC,
                    std::function<double(double,double)> f)
{
#define Lb        BdC.Lboundary
#define Rb        BdC.Rboundary
#define Ub        BdC.Uboundary
#define Db        BdC.Dboundary
#define Lm        BdC.LMix
#define Rm        BdC.RMix
#define Um        BdC.UMix
#define Dm        BdC.DMix
#define Dirichlet boundtype::Dirichlet
#define Neumann   boundtype::Neumann

    int    n = size;
    int    m = n*n;
    double h = 1.0/(n-1);
    std::vector<double> A(m*m); // row major
    double b[m];

    // for the points in the first row
    if(Lm.first==0 && Dm.first==0){
        A[0] = 1;
        b[0] = Db.first(0,0);
    }else{
        A[0] = 4*Lm.first*Dm.first+2*h*Lm.second*Dm.first+2*h*Lm.first*Dm.second;
        A[1] = -2*Lm.first*Dm.first;
        A[n] = -2*Lm.first*Dm.first;
        b[0] = Lm.first*Dm.first*h*h*f(0,0)+2*h*Dm.first*Lb.first(0,0)+2*h*Lm.first*Db.first(0,0);
    }
    if(Rm.first==0 && Dm.first==0){
        A[(n-1)*(m+1)] = 1;
        b[n-1] = Db.first(1,0);
    }else{
        A[(n-1)*(m+1)] = 4*Rm.first*Dm.first+2*h*Rm.second*Dm.first+2*h*Rm.first*Dm.second;
        A[(n-1)*(m+1)-1] = -2*Rm.first*Dm.first;
        A[(n-1)*(m+1)+n] = -2*Rm.first*Dm.first;
        b[n-1] = Rm.first*Dm.first*h*h*f(1,0)+2*h*Dm.first*Rb.first(1,0)+2*h*Rm.first*Db.first(1,0);
    }
    for(int j=1; j<n-1; j++){
        A[j*m + j] = 4*Dm.first+2*h*Dm.second;
        A[j*m+j-1] = -Dm.first;
        A[j*m+j+1] = -Dm.first;
        A[j*m+j+n] = -2*Dm.first;
        b[j] = Dm.first*h*h*f(j*h,0)+2*h*Db.first(j*h,0);
    }
    // for the points in the middlle n-2 rows 
    for(int i=1; i<n-1; i++){
        for(int j=1; j<n-1; j++){
            int index = (i*n+j)*m + i*n+j;
            A[index]   =  4;
            A[index-n] = -1;
            A[index+n] = -1;
            A[index-1] = -1;
            A[index+1] = -1;
            b[i*n+j]   = h*h*f(j*h,i*h);
        }
    }
    for(int i=1; i<n-1; i++){
        int index = i*n*m + i*n;
        A[index] = 4*Lm.first+2*h*Lm.second;
        A[index-n] = -Lm.first;
        A[index+n] = -Lm.first;
        A[index+1] = -2*Lm.first;
        b[i*n] = Lm.first*h*h*f(0,i*h)+2*h*Lb.first(0,i*h);
    }
    for(int i=1; i<n-1; i++){
        int index = (i*n+n-1)*(m+1);
        A[index] = 4*Rm.first+2*h*Rm.second;
        A[index-n] = -Rm.first;
        A[index+n] = -Rm.first;
        A[index-1] = -2*Rm.first;
        b[i*n+n-1] = Rm.first*h*h*f(1,i*h)+2*h*Rb.first(1,i*h);
    }
    // for the points in the last row
    if(Lm.first==0 && Um.first==0){
        A[(n*n-n)*(m+1)] = 1;
        b[n*n-n] = Ub.first(0,1);
    }else{
        A[(n*n-n)*(m+1)] = 4*Lm.first*Um.first+2*h*Lm.second*Um.first+2*h*Lm.first*Um.second;
        A[(n*n-n)*(m+1)-n] = -2*Lm.first*Um.first;
        A[(n*n-n)*(m+1)+1] = -2*Lm.first*Um.first;
        b[n*n-n] = Lm.first*Um.first*h*h*f(0,1)+2*h*Um.first*Lb.first(0,1)+2*h*Lm.first*Ub.first(0,1);
    }
    if(Rm.first==0 && Um.first==0){
        A[(n*n-1)*m+n*n-1] = 1;
        b[n*n-1] = Ub.first(1,1);
    }else{
        A[(n*n-1)*m+n*n-1] = 4*Rm.first*Um.first+2*h*Rm.second*Um.first+2*h*Rm.first*Um.second;
        A[(n*n-1)*m+n*n-1-n] = -2*Rm.first*Um.first;
        A[(n*n-1)*m+n*n-1-1] = -2*Rm.first*Um.first;
        b[n*n-1] = Rm.first*Um.first*h*h*f(1,1)+2*h*Um.first*Rb.first(1,1)+2*h*Rm.first*Ub.first(1,1);
    }
    for(int i=1; i<n-1; i++){
        int index = (n*n-n+i)*m + n*n-n+i;
        A[index] = 4*Um.first+2*h*Um.second;
        A[index-1] = -Um.first;
        A[index+1] = -Um.first;
        A[index-n] = -2*Um.first;
        b[n*n-n+i] = Um.first*h*h*f(i*h,1)+2*h*Ub.first(i*h,1);
    }
    // check if all the four boundary conditions are Neumann
    if(Lm.second==0 && Rm.second==0 && Um.second==0 && Dm.second==0){
        for(int i=1; i<m; i++){
            A[i] = 0;
        }
        A[0] = 1; b[0] = BdC.value_at_00;
    }

    int ipiv[m];
    LAPACKE_dgesv(LAPACK_ROW_MAJOR,m,1,&A[0],m,ipiv,b,1);

    for(int i=0; i<m; i++){
        points[i].value = b[i];
    }

#undef Neumann 
#undef Dirichlet
#undef Dm
#undef Um
#undef Rm
#undef Lm
#undef Db
#undef Ub
#undef Rb
#undef Lb
}

void grid<domin::rmdisk>::BVP_solver(
                    BoundaryCondi<domin::rmdisk> BdC,
                    std::function<double(double,double)> f)
{
#define Lb BdC.Lboundary
#define Rb BdC.Rboundary
#define Ub BdC.Uboundary
#define Db BdC.Dboundary
#define Cb BdC.Cirboundary
#define Lm BdC.LMix
#define Rm BdC.RMix
#define Um BdC.UMix
#define Dm BdC.DMix
#define Cm BdC.CMix
#define Dirichlet boundtype::Dirichlet
#define Neumann boundtype::Neumann

    int    n = points.size();
    double h = 1.0/(size-1);
    std::vector<double> A(n*n);
    double B[n];
    for(int i=0; i<n; i++){
        double x = points[i]._x;
        double y = points[i]._y;
        switch (points[i].pos_type)
        {
        case innerpoint::corner :{
            if(x==0 && y==0){
                if(Lm.first==0 && Dm.first==0){
                    A[i*n + i] = 1;
                    B[i] = Db.first(x,y);
                }else{
                    double a = Lm.first, b = Lm.second;
                    double c = Dm.first, d = Dm.second;
                    auto right = right_pnt(x,y);
                    auto upper = upper_pnt(x,y);
                    double thetaH = points[i].distance(right);
                    double alphaH = points[i].distance(upper);
                    A[i*n+i] = 2*a*c*(pow(thetaH,2)+pow(alphaH,2)) +
                                2*b*c*thetaH*pow(alphaH,2) + 2*a*d*pow(thetaH,2)*alphaH;
                    A[i*n+ret_order(right)] = -2*a*c*pow(alphaH,2);
                    A[i*n+ret_order(upper)] = -2*a*c*pow(thetaH,2);
                    B[i] = a*c*pow(thetaH,2)*pow(alphaH,2)*f(x,y) + 
                            2*c*thetaH*pow(alphaH,2)*Lb.first(x,y) + 
                            2*a*alphaH*pow(thetaH,2)*Db.first(x,y);
                }
            }else if(x==0 && y==1){
                if(Lm.first==0 && Um.first==0){
                    A[i*n+i] = 1;
                    B[i] = Ub.first(x,y);
                }else{
                    double a = Lm.first, b = Lm.second;
                    double c = Um.first, d = Um.second;
                    auto right = right_pnt(x,y);
                    auto down = down_pnt(x,y);
                    double thetaH = points[i].distance(right);
                    double alphaH = points[i].distance(down);
                    A[i*n+i] = 2*a*c*(pow(thetaH,2)+pow(alphaH,2)) +
                                2*b*c*thetaH*pow(alphaH,2) + 2*a*d*pow(thetaH,2)*alphaH;
                    A[i*n+ret_order(right)] = -2*a*c*pow(alphaH,2);
                    A[i*n+ret_order(down)] = -2*a*c*pow(thetaH,2);
                    B[i] = a*c*pow(thetaH,2)*pow(alphaH,2)*f(x,y) + 
                            2*c*thetaH*pow(alphaH,2)*Lb.first(x,y) + 
                            2*a*alphaH*pow(thetaH,2)*Ub.first(x,y);
                }
            }else if(x==1 && y==0){
                if(Rm.first==0 && Dm.first==0){
                    A[i*n+i] = 1;
                    B[i] = Db.first(x,y);
                }else{
                    double a = Rm.first, b = Rm.second;
                    double c = Dm.first, d = Dm.second;
                    auto left = left_pnt(x,y);
                    auto upper = upper_pnt(x,y);
                    double thetaH = points[i].distance(left);
                    double alphaH = points[i].distance(upper);
                    A[i*n+i] = 2*a*c*(pow(thetaH,2)+pow(alphaH,2)) +
                                2*b*c*thetaH*pow(alphaH,2) + 2*a*d*pow(thetaH,2)*alphaH;
                    A[i*n+ret_order(left)] = -2*a*c*pow(alphaH,2);
                    A[i*n+ret_order(upper)] = -2*a*c*pow(thetaH,2);
                    B[i] = a*c*pow(thetaH,2)*pow(alphaH,2)*f(x,y) + 
                            2*c*thetaH*pow(alphaH,2)*Rb.first(x,y) + 
                            2*a*alphaH*pow(thetaH,2)*Db.first(x,y);
                }
            }else{
                if(Rm.first==0 && Um.first==0){
                    A[i*n+i] = 1;
                    B[i] = Ub.first(x,y);
                }else{
                    double a = Rm.first, b = Rm.second;
                    double c = Um.first, d = Um.second;
                    auto left = left_pnt(x,y);
                    auto down = down_pnt(x,y);
                    double thetaH = points[i].distance(left);
                    double alphaH = points[i].distance(down);
                    A[i*n+i] = 2*a*c*(pow(thetaH,2)+pow(alphaH,2)) +
                                2*b*c*thetaH*pow(alphaH,2) + 2*a*d*pow(thetaH,2)*alphaH;
                    A[i*n+ret_order(left)] = -2*a*c*pow(alphaH,2);
                    A[i*n+ret_order(down)] = -2*a*c*pow(thetaH,2);
                    B[i] = a*c*pow(thetaH,2)*pow(alphaH,2)*f(x,y) + 
                            2*c*thetaH*pow(alphaH,2)*Rb.first(x,y) + 
                            2*a*alphaH*pow(thetaH,2)*Ub.first(x,y);
                }
            }
            break;}
        case innerpoint::insquare :{
            auto left  = left_pnt(x,y); auto right = right_pnt(x,y);
            auto down  = down_pnt(x,y); auto upper = upper_pnt(x,y);
            double lef_dis = points[i].distance(left);
            double rig_dis = points[i].distance(right);
            double upp_dis = points[i].distance(upper);
            double dow_dis = points[i].distance(down);
            A[i*n + i] = (lef_dis+rig_dis)/(lef_dis*rig_dis*(lef_dis+rig_dis)) + 
                            (upp_dis+dow_dis)/(upp_dis*dow_dis*(upp_dis+dow_dis));
            A[i*n + ret_order(left)]  = -rig_dis/(lef_dis*rig_dis*(lef_dis+rig_dis));
            A[i*n + ret_order(right)] = -lef_dis/(lef_dis*rig_dis*(lef_dis+rig_dis));
            A[i*n + ret_order(upper)] = -dow_dis/(upp_dis*dow_dis*(upp_dis+dow_dis));
            A[i*n + ret_order(down)]  = -upp_dis/(upp_dis*dow_dis*(upp_dis+dow_dis));
            B[i] = 0.5*f(x,y);
            break;}
        case innerpoint::leftbound :{
            double a = Lm.first, b = Lm.second;
            auto right = right_pnt(x,y);
            auto upper = upper_pnt(x,y); auto down = down_pnt(x,y);
            double rig_dis = points[i].distance(right);
            double upp_dis = points[i].distance(upper);
            double dow_dis = points[i].distance(down);
            A[i*n+i] = 1.0*a/pow(rig_dis,2) + a*1.0/(upp_dis*dow_dis)
                        + b*1.0/rig_dis;
            A[i*n+ret_order(right)] = -a*1.0/pow(rig_dis,2);
            A[i*n+ret_order(upper)] = -a*1.0/(upp_dis*(upp_dis+dow_dis));
            A[i*n+ret_order(down)]  = -a*1.0/(dow_dis*(upp_dis+dow_dis));
            B[i] = 0.5*a*f(x,y) + Lb.first(x,y)*1.0/rig_dis;
            break;}
        case innerpoint::rightbound :{
            double a = Rm.first, b = Rm.second;
            auto left = left_pnt(x,y);
            auto upper = upper_pnt(x,y); auto down = down_pnt(x,y);
            double lef_dis = points[i].distance(left);
            double upp_dis = points[i].distance(upper);
            double dow_dis = points[i].distance(down);
            A[i*n + i] = a*1.0/pow(lef_dis,2) + a*1.0/(upp_dis*dow_dis)
                            + b*1.0/lef_dis;
            A[i*n + ret_order(left)]  = -a*1.0/pow(lef_dis,2);
            A[i*n + ret_order(upper)] = -a*1.0/(upp_dis*(upp_dis+dow_dis));
            A[i*n + ret_order(down)]  = -a*1.0/(dow_dis*(upp_dis+dow_dis));
            B[i] = 0.5*a*f(x,y) + Rb.first(x,y)*1.0/lef_dis;
            break;}
        case innerpoint::upbound :{
            double a = Um.first, b = Um.second;
            auto down = down_pnt(x,y);
            auto left = left_pnt(x,y); auto right = right_pnt(x,y);
            double dow_dis = points[i].distance(down);
            double lef_dis = points[i].distance(left);
            double rig_dis = points[i].distance(right);
            A[i*n + i] = a*1.0/pow(dow_dis,2) + a*1.0/(lef_dis*rig_dis)
                            + b*1.0/dow_dis;
            A[i*n + ret_order(down)]  = -a*1.0/pow(dow_dis,2);
            A[i*n + ret_order(left)]  = -a*1.0/(lef_dis*(lef_dis+rig_dis));
            A[i*n + ret_order(right)] = -a*1.0/(rig_dis*(lef_dis+rig_dis));
            B[i] = 0.5*a*f(x,y) + Ub.first(x,y)*1.0/dow_dis;
            break;}
        case innerpoint::downbound :{
            double a = Dm.first, b = Dm.second;
            auto upper = upper_pnt(x,y);
            auto left = left_pnt(x,y); auto right = right_pnt(x,y);
            double upp_dis = points[i].distance(upper);
            double lef_dis = points[i].distance(left);
            double rig_dis = points[i].distance(right);
            A[i*n + i] = a*1.0/pow(upp_dis,2) + a*1.0/(lef_dis*rig_dis)
                            + b*1.0/upp_dis;
            A[i*n + ret_order(upper)] = -a*1.0/pow(upp_dis,2);
            A[i*n + ret_order(left)]  = -a*1.0/(lef_dis*(lef_dis+rig_dis));
            A[i*n + ret_order(right)] = -a*1.0/(rig_dis*(lef_dis+rig_dis));
            B[i] = 0.5*a*f(x,y) + Db.first(x,y)*1.0/upp_dis;
            break;}
        case innerpoint::ondisk :{
            if(Cb.second == Dirichlet){
                A[i*n + i] = 1;
                B[i] = Cb.first(x,y);
            }else{
                auto insection = intersection_with_grid(D._x0,D._y0,x,y);
            #define OUT_OF_BOUND                                \
                double alpha = insection.distance(points[i]);   \
                double beta  = insection.distance(dot1);        \
                double sigma = insection.distance(dot2);        \
                A[i*n + i]               = -1;                  \
                A[i*n + ret_order(dot1)] = -sigma/(beta-sigma); \
                A[i*n + ret_order(dot2)] = beta/(beta-sigma);   \
                B[i] = alpha*Cb.first(x,y);
            #define INSIDE_BOUND                                \
                double alpha = insection.distance(points[i]);   \
                double beta  = insection.distance(dot1);        \
                double sigma = insection.distance(dot2);        \
                A[i*n + i]               = -1;                  \
                A[i*n + ret_order(dot1)] = sigma/(beta+sigma);  \
                A[i*n + ret_order(dot2)] = beta/(beta+sigma);   \
                B[i] = alpha*Cb.first(x,y);
                switch(insection.pos_type)
                {
                case innerpoint::onXaxis :
                    if(insection._x < 0){
                        auto dot1 = find(0,insection._y);
                        auto dot2 = right_pnt(dot1._x,dot1._y);
                        OUT_OF_BOUND
                    }else if(insection._x > 1){
                        auto dot1 = find(1,insection._y);
                        auto dot2 = left_pnt(dot1._x,dot1._y);
                        OUT_OF_BOUND
                    }else{
                        int axis1 = static_cast<int>(insection._x/h);
                        int axis2 = axis1+1;
                        auto dot1 = find(axis1*h, insection._y);
                        auto dot2 = find(axis2*h, insection._y);
                        INSIDE_BOUND
                    }
                    break;
                case innerpoint::onYaxis :
                    if(insection._y < 0){
                        auto dot1 = find(insection._x,0);
                        auto dot2 = upper_pnt(dot1._x,dot1._y);
                        OUT_OF_BOUND
                    }else if(insection._y > 1){
                        auto dot1 = find(insection._x,1);
                        auto dot2 = down_pnt(dot1._x,dot1._y);
                        OUT_OF_BOUND
                    }else{
                        int axis1 = static_cast<int>(insection._y/h);
                        int axis2 = axis1+1;
                        auto dot1 = find(insection._x, axis1*h);
                        auto dot2 = find(insection._x, axis2*h);
                        INSIDE_BOUND
                    }
                    break;
                }
            #undef INSIDE_BOUND
            #undef OUT_OF_BOUND
            }
            break;}
        }
    }

    int ipiv[n];
    LAPACKE_dgesv(LAPACK_ROW_MAJOR,n,1,&A[0],n,ipiv,B,1);

    for(int i=0; i<n; i++){
        points[i].value = B[i];
    }

#undef Neumann
#undef Dirichlet
#undef Cm
#undef Dm
#undef Um
#undef Rm
#undef Lm
#undef Cb
#undef Db
#undef Ub
#undef Rb
#undef Lb
}

#endif