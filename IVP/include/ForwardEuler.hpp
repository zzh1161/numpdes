#ifndef __FORWARD_EULER_HPP__
#define __FORWARD_EULER_HPP__

#include"TimeIntegrator.hpp"

template<int dim>
class ForwardEulerSolver: public TimeIntegrator<dim>
{
public:
    using Func = std::function<Vec<real,dim>(Vec<real,dim>,real)>;

    ForwardEulerSolver(int steps): steps_(steps){}
    ForwardEulerSolver(const ForwardEulerSolver &) = delete;
    ForwardEulerSolver& operator=(const ForwardEulerSolver &) = delete;
    ~ForwardEulerSolver() = default;

    Vec<real,dim> solve(const Vec<real,dim> &U0, Func f, real T, int p_=0) const
    {
        #define U TimeIntegrator<dim>::result_
        std::vector<Vec<real,dim>>().swap(U);
        real k = 1.0*T/steps_;
        U.push_back(U0);
        for(int n=0; n<steps_; n++){
            U.push_back(U[n] + f(U[n],n*k)*k);
        }
        #undef U
        return *(TimeIntegrator<dim>::result_.end()-1);
    }

    real solution_error(const Vec<real,dim> &U0, Func f, real T, int p_=0)
    {
        #define U TimeIntegrator<dim>::result_
        Vec<real,dim> res1;
        if(U.size() == 0)
            res1 = solve(U0, f, T, p_);
        else
            res1 = *(U.end()-1);
        std::vector<Vec<real,dim>>().swap(U);
        real k = 1.0*T/steps_;
        U.push_back(U0);
        for(int n=0; n<steps_; n++){
            auto u1 = U[n] + f(U[n],n*k)*k;
            auto u_half = U[n] + f(U[n], n*k)*(0.5*k);
            auto u2 = u_half + f(u_half, n*k+0.5*k)*(0.5*k);
            U.push_back(u2*2-u1);
        }
        auto res2 = *(U.end()-1);
        real error = norm(res1-res2, 0);
        #undef U
        return error;
    }

private:
    int steps_;
};

template<int dim>
std::unique_ptr<TimeIntegrator<dim>> createForwardEuler(int steps=100)
{
    return std::make_unique<ForwardEulerSolver<dim>>(steps);
}


#endif