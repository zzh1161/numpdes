#ifndef __DP54EMBED_RK_HPP__
#define __DP54EMBED_RK_HPP__

#include"TimeIntegrator.hpp"

template<int dim>
class DP54embedRKsolver : public TimeIntegrator<dim>
{
public:
    using Func = std::function<Vec<real,dim>(Vec<real,dim>,real)>;

    DP54embedRKsolver(int steps): steps_(steps){}
    DP54embedRKsolver(const DP54embedRKsolver &) = delete;
    DP54embedRKsolver& operator=(const DP54embedRKsolver &) = delete;
    ~DP54embedRKsolver() = default;

    Vec<real,dim> solve(const Vec<real,dim> &U0, Func f, real T, int p_=0) const
    {
        #define U TimeIntegrator<dim>::result_
        std::vector<Vec<real,dim>>().swap(U);
        real k = 1.0*T/steps_;
        U.push_back(U0);
        real t = 0; int n = 0; real epsilon = 1e-6;
        // for(int n=0; n<steps_; n++){
        //     auto y1 = f(U[n], n*k);
        //     auto y2 = f(U[n]+y1*(0.2*k), n*k+0.2*k);
        //     auto y3 = f(U[n]+y1*(3*k/40)+y2*(9*k/40), n*k+0.3*k);
        //     auto y4 = f(U[n]+y1*(44*k/45)-y2*(56*k/15)+y3*(32*k/9), n*k+0.8*k);
        //     auto y5 = f(U[n] + y1*(19372*k/6561) - y2*(25360.0*k/2187)
        //                 + y3*(64448*k/6561) - y4*(212*k/729), n*k+(8*k/9));
        //     auto y6 = f(U[n] + y1*(9017*k/3168) - y2*(355*k/33) + y3*(46732*k/5247)
        //                 + y4*(49*k/176) - y5*(5103*k/18656), n*k+k);
        //     auto y7 = f(U[n] + y1*(35*k/384) + y3*(500.0*k/1113) + y4*(125*k/192)
        //                 - y5*(2187*k/6784) + y6*(11*k/84), n*k+k);
        //     U.push_back(U[n] + y1*(35*k/384) + y3*(500.0*k/1113) + y4*(125*k/192)
        //                 - y5*(2187*k/6784) + y6*(11*k/84));
        // }
        while(t < T){
            k = std::min(k, T-t);
            auto y1 = f(U[n], t);
            auto y2 = f(U[n]+y1*(0.2*k), t+0.2*k);
            auto y3 = f(U[n]+y1*(3*k/40)+y2*(9*k/40), t+0.3*k);
            auto y4 = f(U[n]+y1*(44*k/45)-y2*(56*k/15)+y3*(32*k/9), t+0.8*k);
            auto y5 = f(U[n] + y1*(19372*k/6561) - y2*(25360.0*k/2187)
                        + y3*(64448*k/6561) - y4*(212*k/729), t+(8*k/9));
            auto y6 = f(U[n] + y1*(9017*k/3168) - y2*(355*k/33) + y3*(46732*k/5247)
                        + y4*(49*k/176) - y5*(5103*k/18656), t+k);
            auto y7 = f(U[n] + y1*(35*k/384) + y3*(500.0*k/1113) + y4*(125*k/192)
                        - y5*(2187*k/6784) + y6*(11*k/84), t+k);
            auto u1 = U[n] + y1*(35*k/384) + y3*(500.0*k/1113) + y4*(125*k/192)
                      - y5*(2187*k/6784) + y6*(11*k/84);
            auto u2 = U[n] + y1*(5179*k/57600) + y3*(7571*k/16695) + y4*(393*k/640)
                      - y5*(92097*k/339200) + y6*(187*k/2100) + y7*(k/40);
            real R = norm(u1-u2,0)/k;
            real delta = 0.84*pow(epsilon/R, 0.25);
            if(R <= epsilon){
                t += k;
                U.push_back(u1);
                n++;
            }
            k = k*delta;
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
std::unique_ptr<TimeIntegrator<dim>> createDP54RK(int steps=100)
{
    return std::make_unique<DP54embedRKsolver<dim>>(steps);
}


#endif