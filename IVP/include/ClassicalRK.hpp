#ifndef __CLASSICAL_RK_HPP__
#define __CLASSICAL_RK_HPP__

#include"TimeIntegrator.hpp"

template<int dim>
class ClassicalRKsolver : public TimeIntegrator<dim>
{
public:
    using Func = std::function<Vec<real,dim>(Vec<real,dim>,real)>;

    ClassicalRKsolver(int steps): steps_(steps){}
    ClassicalRKsolver(const ClassicalRKsolver &) = delete;
    ClassicalRKsolver& operator=(const ClassicalRKsolver &) = delete;
    ~ClassicalRKsolver() = default;

    Vec<real,dim> solve(const Vec<real,dim> &U0, Func f, real T, int p_=0) const
    {
        #define U TimeIntegrator<dim>::result_
        std::vector<Vec<real,dim>>().swap(U);
        real k = 1.0*T/steps_;
        U.push_back(U0);
        for(int n=0; n<steps_; n++){
            auto y1 = f(U[n], n*k);
            auto y2 = f(U[n]+y1*(0.5*k), n*k+0.5*k);
            auto y3 = f(U[n]+y2*(0.5*k), n*k+0.5*k);
            auto y4 = f(U[n]+y3*k, n*k+k);
            U.push_back(U[n] + (y1+y2*2+y3*2+y4)*(1.0*k/6));
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
std::unique_ptr<TimeIntegrator<dim>> createClassiRK(int steps=100)
{
    return std::make_unique<ClassicalRKsolver<dim>>(steps);
}


#endif