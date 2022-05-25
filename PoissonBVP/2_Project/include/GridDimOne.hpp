/**************************************************************
 > @Description    : MultiGrid methods in one-dimensional grid
 > @Version        : 1.0
 > @Author         : zhang-zh
 > @Date           : 2022-04-12 22:29
 > @LastEditTime   : 2022-04-26 21:29
**************************************************************/
#ifndef __PRO2_GRID_DIMONE_HPP__
#define __PRO2_GRID_DIMONE_HPP__

#include"Grid.hpp"

template<>
class Grid<1,Domin::regular>
{
public:
    explicit Grid(const int n, BoundaryCondition<1> &B, 
                        std::function<real(real)> f){
        assert((n > 0) && (n&(n-1))==0 && n>=MINGrid);
        real h = 1.0/n;
        size_ = n;
        BoundCond_ = B;
        points_.resize(n+1);
        RightHandSide_.resize(n+1);
        RightHandFunc_ = f;
        for(int i=1; i<n; i++){
            // (-1 2 -1)*1/h^2
            RightHandSide_[i] = f(i*h);
        }
        // (3a/2h +b)U_0 - 2a/h U_1 + a/2h U_2 = sigma
        RightHandSide_[0] = BoundCond_.leftBound.first; 
        // (3c/2h +d)U_n - 2c/h U_(n-1) + c/2h U_(n-2) = sigma
        RightHandSide_[n] = BoundCond_.rightBound.first;
    }

private:
    /**
     * @param  w       weighted Jacobi's coefficient
     * @param  Rest    restriction operator
     * @param  Inter   interpolation operator
     */
    std::vector<real> V_cycle(real w, std::vector<real> v, std::vector<real> f, int nu1,
            int nu2, Operators::Restriction Rest, Operators::Interpolation Inter);
    std::vector<real> FMG(real w, std::vector<real> f, int nu1, 
            int nu2, Operators::Restriction Rest, Operators::Interpolation Inter);

public:
    void Poisson_BVP_Multigrid(MultigridMethod method, real w = 1, int nu1 = 3, int nu2 = 3, 
                               Operators::Restriction Rest = Operators::fullWeighting,
                               Operators::Interpolation Inter = Operators::linear,
                               int max_iter = 100, real epsilon = 0.0001);

    real analytic_error(std::function<real(real)> f) const{
        std::vector<real> err(size_+1);
        for(int i=0; i<=size_; i++)
            err[i] = abs(points_[i]-f(i*1.0/size_));
        auto Max = std::max_element(err.begin(), err.end());
        return *Max;
    }

    real discrete_error() const{
        auto n = size_;
        std::vector<real> u = RightHandSide_;
        std::vector<real> e(n+1);
        std::vector<real> A((n+1)*(n+1));
        auto a = BoundCond_.leftMixedCoeffi.first; 
        auto b = BoundCond_.leftMixedCoeffi.second;
        auto c = BoundCond_.rightMixedCoeffi.first;
        auto d = BoundCond_.rightMixedCoeffi.second;
        A[0] = 1.5*a*n+b; A[1] = -2*a*n; A[2] = 0.5*a*n;
        A[(n+1)*(n+1)-1] = 1.5*c*n+d;
        A[(n+1)*(n+1)-2] = -2*c*n;
        A[(n+1)*(n+1)-3] = 0.5*c*n;
        for(int i=1; i<n; i++){
            A[i*(n+1) + i]   = 2*n*n;
            A[i*(n+1) + i-1] = -1*n*n;
            A[i*(n+1) + i+1] = -1*n*n;
        }
        int ipiv[n+1];
        LAPACKE_dgesv(LAPACK_ROW_MAJOR, n+1, 1, &A[0], n+1, ipiv, &u[0], 1);
        for(int i=0; i<=n; i++)
            e[i] = abs(points_[i]-u[i]);
        // real norm_e = sqrt(std::accumulate(e.begin(),
        //                                    e.end(),
        //                                    0.0,
        //                                    [](real x, real y){return x+y*y;}));
        // real norm_u = sqrt(std::accumulate(u.begin(),
        //                                    u.end(),
        //                                    0.0,
        //                                    [](real x, real y){return x+y*y;}));
        // return norm_e/norm_u;
        return *std::max_element(e.begin(), e.end());
    }

    void error_convergence(std::function<real(real)> f);

    void display() const{
        real h = 1.0/size_;
        std::cout << "x = [";
        for(int i=0; i<=size_; i++)
            std::cout << i*h << ", ";
        std::cout << "];" << std::endl;
        std::cout << "y = [";
        for(int i=0; i<=size_; i++)
            std::cout << points_[i] << ", ";
        std::cout << "];" << std::endl;
    }

private:
    int size_;
    std::vector<real>         points_;
    BoundaryCondition<1>      BoundCond_;
    std::function<real(real)> RightHandFunc_;
    std::vector<real>         RightHandSide_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////
/************************************** V-cycle & FMG ********************************************/
///////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<real> Grid<1,Domin::regular>::V_cycle(real w, std::vector<real> v, std::vector<real> f, 
                    int nu1, int nu2, Operators::Restriction Rest, Operators::Interpolation Inter){
    int n = v.size()-1;
    assert((n > 0) && (n&(n-1))==0 && n>=MINGrid);
    assert(v.size()==n+1 && f.size()==n+1);
    real h = 1.0/n;
    auto a = BoundCond_.leftMixedCoeffi.first;
    auto b = BoundCond_.leftMixedCoeffi.second;
    auto c = BoundCond_.rightMixedCoeffi.first;
    auto d = BoundCond_.rightMixedCoeffi.second;
    // Relax nu1 times
    std::vector<real> temp(n+1);
    for(int times=0; times<nu1; times++){
        for(int i=1; i<n; i++){
            temp[i] = 0.5*w*v[i-1] + 0.5*w*v[i+1] + (1-w)*v[i] + 0.5*w*h*h*f[i];
        }
        temp[0] = (1-w)*v[0] + (4*a*n*w/(3*a*n+2*b))*v[1] -
                  (a*n*w/(3*a*n+2*b))*v[2] + (2*w/(3*a*n+2*b))*f[0];
        temp[n] = (1-w)*v[n] + (4*c*n*w/(3*c*n+2*d))*v[n-1] -
                  (c*n*w/(3*c*n+2*d))*v[n-2] + (2*w/(3*c*n+2*d))*f[n];
        v = temp;
    }
    // Recursion
    if(n != MINGrid){
        // Restriction
        std::vector<real> bak_f(n+1);
        std::vector<real> f_2(n/2+1);
        std::vector<real> v_2(n/2+1); //v_2[0]=RightHandSide_[0]; v_2[n/2]=RightHandSide_[size_];
        bak_f[0] = f[0] - (1.5*a*n+b)*v[0]+2*a*n*v[1]-0.5*a*n*v[2];
        bak_f[n] = f[n] - (1.5*c*n+d)*v[n]+2*c*n*v[n-1]-0.5*c*n*v[n-2];
        for(int i=1; i<n; i++){
            bak_f[i] = f[i] + n*n*v[i-1]-2*n*n*v[i]+n*n*v[i+1];
        }
        if(Rest == Operators::fullWeighting){
            f_2[0] = bak_f[0]; f_2[n/2] = bak_f[n];
            // f_2[0]   = 0.5*bak_f[0]+0.25*bak_f[1];
            // f_2[n/2] = 0.5*bak_f[n]+0.25*bak_f[n-1];
            for(int i=1; i<n/2; i++)
                f_2[i] = 0.25*bak_f[2*i-1]+0.5*bak_f[2*i]+0.25*bak_f[2*i+1];
        }else if(Rest == Operators::injection){
            for(int i=0; i<=n/2; i++)
                f_2[i] = bak_f[2*i];
        }
        // Recursion
        v_2 = V_cycle(w, v_2, f_2, nu1, nu2, Rest, Inter);
        // Prolongation
        if(Inter == Operators::linear){
            for(int i=0; i<=n; i+=2)
                v[i] += v_2[i/2];
            // for(int i=3; i<n-1; i+=2)
            //     v[i] += 0.5*v_2[(i-1)/2] + 0.5*v_2[(i+1)/2];
            // v[1] += 0.5*v_2[1]; v[n-1] += 0.5*v_2[n/2-1];
            for(int i=1; i<n; i+=2)
                v[i] += 0.5*v_2[(i-1)/2] + 0.5*v_2[(i+1)/2];
        }else if(Inter == Operators::quadratic){
            for(int i=0; i<=n; i+=2)
                v[i] += v_2[i/2];
            for(int i=2; i<=n-2; i+=4){
                v[i-1] += v_2[i/2-1]*0.375 + v_2[i/2]*0.75 - v_2[i/2+1]*0.125;
                v[i+1] += -v_2[i/2-1]*0.125 + v_2[i/2]*0.75 + v_2[i/2+1]*0.375;
            }
        }
    }
    // Relax nu2 times
    for(int times=0; times<nu2; times++){
        for(int i=1; i<n; i++){
            temp[i] = 0.5*w*v[i-1] + 0.5*w*v[i+1] + (1-w)*v[i] + 0.5*w*h*h*f[i];
        }
        temp[0] = (1-w)*v[0] + (4*a*n*w/(3*a*n+2*b))*v[1] -
                  (a*n*w/(3*a*n+2*b))*v[2] + (2*w/(3*a*n+2*b))*f[0];
        temp[n] = (1-w)*v[n] + (4*c*n*w/(3*c*n+2*d))*v[n-1] -
                  (c*n*w/(3*c*n+2*d))*v[n-2] + (2*w/(3*c*n+2*d))*f[n];
        v = temp;
    }

    return v;
}

std::vector<real> Grid<1,Domin::regular>::FMG(real w, std::vector<real> f, 
            int nu1, int nu2, Operators::Restriction Rest, Operators::Interpolation Inter){   
    int n  = f.size()-1;
    real h = 1.0/n;
    assert((n > 0) && (n&(n-1))==0 && n>=MINGrid);
    // Judgement
    if(n == MINGrid){
        std::vector<real> v(n+1);
        return V_cycle(w,v,f,nu1,nu2,Rest,Inter);
    }
    // Restriction and Recursion
    std::vector<real> f_2(n/2+1);
    if(Rest == Operators::fullWeighting){
        f_2[0] = f[0]; f_2[n/2] = f[n];
        // f_2[0]   = 0.5*f[0]+0.25*f[1];
        // f_2[n/2] = 0.5*f[n]+0.25*f[n-1];
        for(int i=1; i<n/2; i++)
            f_2[i] = 0.25*f[2*i-1]+0.5*f[2*i]+0.25*f[2*i+1];
    }else if(Rest == Operators::injection){
        for(int i=0; i<=n/2; i++)
            f_2[i] = f[2*i];
    }
    auto v_2 = FMG(w,f_2,nu1,nu2,Rest,Inter);
    // prolongation
    std::vector<real> v(n+1);
    if(Inter == Operators::linear){
        for(int i=0; i<=n; i+=2)
            v[i] = v_2[i/2];
        // for(int i=3; i<n-1; i+=2)
        //     v[i] = 0.5*v_2[(i-1)/2] + 0.5*v_2[(i+1)/2];
        // v[1] = 0.5*v_2[1]; v[n-1] = 0.5*v_2[n/2-1];
        for(int i=1; i<n; i+=2)
            v[i] = 0.5*v_2[(i-1)/2] + 0.5*v_2[(i+1)/2];
    }else if(Inter == Operators::quadratic){
        for(int i=0; i<=n; i+=2)
            v[i] = v_2[i/2];
        for(int i=2; i<=n-2; i+=4){
            v[i-1] = v_2[i/2-1]*0.375 + v_2[i/2]*0.75 - v_2[i/2+1]*0.125;
            v[i+1] = -v_2[i/2-1]*0.125 + v_2[i/2]*0.75 + v_2[i/2+1]*0.375;
        }
    }
    // V_cycle
    return V_cycle(w,v,f,nu1,nu2,Rest,Inter);
}

void Grid<1,Domin::regular>::Poisson_BVP_Multigrid(MultigridMethod method,
                                            real w, int nu1, int nu2, 
                                            Operators::Restriction Rest,
                                            Operators::Interpolation Inter,
                                            int max_iter, real epsilon){
    int n = size_;
    bool flag = false;
    real error;
    std::vector<real> v(size_+1);
    std::vector<real> err(size_+1);
    auto a = BoundCond_.leftMixedCoeffi.first; 
    auto b = BoundCond_.leftMixedCoeffi.second;
    auto c = BoundCond_.rightMixedCoeffi.first;
    auto d = BoundCond_.rightMixedCoeffi.second;
    // v[0] = RightHandSide_[0]; v[size_] = RightHandSide_[size_];
    for(int i=0; i<max_iter; i++){
        if(method == MultigridMethod::Vcycle){
            v = std::move(V_cycle(w, std::move(v), RightHandSide_, nu1, nu2, Rest, Inter));
        }else if(method == MultigridMethod::FMG){
            v = std::move(FMG(w, RightHandSide_, nu1, nu2, Rest, Inter));
            flag = true;
        }
        // residual
        err[0] = abs(RightHandSide_[0] - (1.5*a*n+b)*v[0]+2*a*n*v[1]-0.5*a*n*v[2]);
        err[n] = abs(RightHandSide_[n] - (1.5*c*n+d)*v[n]+2*c*n*v[n-1]-0.5*c*n*v[n-2]);
        for(int i=1; i<n; i++){
            err[i] = abs(RightHandSide_[i] + n*n*v[i-1]-2*n*n*v[i]+n*n*v[i+1]);
        }
        // inf-norm
        error = *std::max_element(err.begin(),err.end());
        if(method == MultigridMethod::Vcycle)
            std::cout << "V-cycle执行第" << i+1 << "次，残差为" << error << std::endl;
        else if(method == MultigridMethod::FMG)
            std::cout << "FMG执行1次，残差为" << error << std::endl;
        if(error < epsilon || flag == true) break;
    }

    points_ = std::move(v);
}

void Grid<1>::error_convergence(std::function<real(real)> f){
    int n = size_;
    int max_iter = 100; real epsilon = 1e-8;
    bool flag = false;
    real error;
    std::vector<real> v(size_+1);
    std::vector<real> err(size_+1);
    std::vector<real> solu(size_+1);
    for(int i=0; i<=size_; i++){
        solu[i] = f(i*1.0/size_);
    }
    auto a = BoundCond_.leftMixedCoeffi.first; 
    auto b = BoundCond_.leftMixedCoeffi.second;
    auto c = BoundCond_.rightMixedCoeffi.first;
    auto d = BoundCond_.rightMixedCoeffi.second;
    for(int i=0; i<max_iter; i++){

        v = std::move(V_cycle(2.0/3, std::move(v), RightHandSide_, 3, 3,
                            Operators::fullWeighting, Operators::linear));
            
        // residual
        err[0] = abs(RightHandSide_[0] - (1.5*a*n+b)*v[0]+2*a*n*v[1]-0.5*a*n*v[2]);
        err[n] = abs(RightHandSide_[n] - (1.5*c*n+d)*v[n]+2*c*n*v[n-1]-0.5*c*n*v[n-2]);
        for(int i=1; i<n; i++){
            err[i] = abs(RightHandSide_[i] + n*n*v[i-1]-2*n*n*v[i]+n*n*v[i+1]);
        }
        // inf-norm
        error = *std::max_element(err.begin(),err.end());
        
        std::cout << "V-cycle执行第" << i+1 << "次，残差为" << error;

        for(int i=0; i<=size_; i++)
            err[i] = abs(v[i]-solu[i]);
        std::cout << "; 与解析解的误差为" << *std::max_element(err.begin(),err.end())
                  << std::endl;

        if(error < epsilon) break;
    }
}

#endif