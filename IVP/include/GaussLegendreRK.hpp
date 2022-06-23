#ifndef __GAUSS_LEGENDRE_RK_HPP__
#define __GAUSS_LEGENDRE_RK_HPP__

#include"TimeIntegrator.hpp"

template<int dim>
class GaussLegendreRKsolver : public TimeIntegrator<dim>
{
public:
    using Func = std::function<Vec<real,dim>(Vec<real,dim>,real)>;

    GaussLegendreRKsolver(int steps): steps_(steps){}
    GaussLegendreRKsolver(const GaussLegendreRKsolver &) = delete;
    GaussLegendreRKsolver& operator=(const GaussLegendreRKsolver &) = delete;
    ~GaussLegendreRKsolver() = default;

    Vec<real,dim> solve(const Vec<real,dim> &U0, Func f, real T, int p_=0) const
    {
        #define U TimeIntegrator<dim>::result_
        std::vector<Vec<real,dim>>().swap(U);
        real k = 1.0*T/steps_;
        auto epsilon = 0.0000001;
        U.push_back(U0);
        if(p_ == 1){
            for(int n=0; n<steps_; n++){
                auto y1 = U[n];
                while(true){
                    auto temp = f(U[n]+y1*(0.5*k), n*k+0.5*k);
                    if(norm(y1-temp, 0) < epsilon) break;
                    else y1 = temp;
                }
                U.push_back(U[n]+y1*k);
            }
        }else if(p_ == 2){
            for(int n=0; n<steps_; n++){
                auto y1 = U[n]; auto y2 = U[n];
                while(true){
                    auto temp1 = f(U[n]+y1*(0.25*k)+y2*((3-2*sqrt(3))*k/12),
                                   n*k+((3-sqrt(3))*k/6));
                    auto temp2 = f(U[n]+y1*((3+2*sqrt(3))*k/12)+y2*(0.25*k),
                                   n*k+((3+sqrt(3))*k/6));
                    if(norm(y1-temp1,0)+norm(y2-temp2,0) < epsilon) break;
                    y1 = temp1; y2 = temp2;
                }
                U.push_back(U[n]+y1*(0.5*k)+y2*(0.5*k));
            }
        }else if(p_ == 3){
            for(int n=0; n<steps_; n++){
                auto y1=U[n]; auto y2=U[n]; auto y3=U[n];
                while(true){
                    auto temp1 = f(U[n] + y1*(5*k/36) + y2*((2.0/9-sqrt(15)/15)*k)
                                   + y3*((5.0/36-sqrt(15)/30)*k), n*k+0.1*(5-sqrt(15))*k);
                    auto temp2 = f(U[n] + y1*((5.0/36+sqrt(15)/24)*k) + y2*(2*k/9)
                                   + y3*((5.0/36-sqrt(15)/24)*k), n*k+0.5*k);
                    auto temp3 = f(U[n] + y1*((5.0/36+sqrt(15)/30)*k) + y2*((2.0/9+sqrt(15)/15)*k)
                                   + y3*(5*k/36), n*k+0.1*(5+sqrt(15))*k);
                    if(norm(y1-temp1,0)+norm(y2-temp2,0)+norm(y3-temp3,0)<epsilon) break;
                    y1 = temp1; y2 = temp2; y3 = temp3;
                }
                U.push_back(U[n]+y1*(5*k/18)+y2*(4*k/9)+y3*(5*k/18));
            }
        }else{
            throw std::runtime_error("No such s for Gauss-Legendre method!");
        }
        #undef U
        return *(TimeIntegrator<dim>::result_.end()-1);
    }

    real solution_error(const Vec<real,dim> &U0, Func f, real T, int p_=0)
    {
        int step_temp = steps_;
        Vec<real,dim> res1;
        if(TimeIntegrator<dim>::result_.size() == 0)
            res1 = solve(U0, f, T, p_);
        else 
            res1 = *(TimeIntegrator<dim>::result_.end()-1);
        steps_ = 2*steps_;
        auto res2 = solve(U0, f, T, p_);
        steps_ = step_temp;
        real error = norm(res1-res2, 0);
        return error;
    }

private:
    int steps_;
};

template<int dim>
std::unique_ptr<TimeIntegrator<dim>> createGaussLege(int steps=100)
{
    return std::make_unique<GaussLegendreRKsolver<dim>>(steps);
}


#endif