#ifndef __F45EMBED_RK_HPP__
#define __F45EMBED_RK_HPP__

#include"TimeIntegrator.hpp"

template<int dim>
class F45embedRKsolver : public TimeIntegrator<dim>
{
public:
    using Func = std::function<Vec<real,dim>(Vec<real,dim>,real)>;

    F45embedRKsolver(int steps): steps_(steps){}
    F45embedRKsolver(const F45embedRKsolver &) = delete;
    F45embedRKsolver& operator=(const F45embedRKsolver &) = delete;
    ~F45embedRKsolver() = default;

    Vec<real,dim> solve(const Vec<real,dim> &U0, Func f, real T, int p_=0) const
    {
        #define U TimeIntegrator<dim>::result_
        std::vector<Vec<real,dim>>().swap(U);
        real k = 1.0*T/steps_;
        U.push_back(U0);
        real t = 0; int n = 0; real epsilon = 1e-6;
        // for(int n=0; n<steps_; n++){
        //     auto y1 = f(U[n], n*k);
        //     auto y2 = f(U[n]+y1*(0.25*k), n*k+0.25*k);
        //     auto y3 = f(U[n]+y1*(3*k/32)+y2*(9*k/32), n*k+(3*k/8));
        //     auto y4 = f(U[n]+y1*(1932*k/2197)-y2*(7200.0*k/2197)+y3*(7296*k/2197), n*k+(12*k/13));
        //     auto y5 = f(U[n] + y1*(439*k/216) - y2*(8*k) + y3*(3680.0*k/513)
        //                 - y4*(845*k/4104), n*k+k);
        //     auto y6 = f(U[n] - y1*(8*k/27) + y2*(2*k) - y3*(3544*k/2565)
        //                 + y4*(1859*k/4104) - y5*(11*k/40), n*k+0.5*k);
        //     U.push_back(U[n] + y1*(25*k/216) + y3*(1408*k/2565)
        //                 + y4*(2197*k/4104) - y5*(0.2*k));
        // }
        while(t < T){
            k = std::min(k, T-t);
            auto y1 = f(U[n], t);
            auto y2 = f(U[n]+y1*(0.25*k), t+0.25*k);
            auto y3 = f(U[n]+y1*(3*k/32)+y2*(9*k/32), t+(3*k/8));
            auto y4 = f(U[n]+y1*(1932*k/2197)-y2*(7200.0*k/2197)+y3*(7296*k/2197), t+(12*k/13));
            auto y5 = f(U[n] + y1*(439*k/216) - y2*(8*k) + y3*(3680.0*k/513)
                        - y4*(845*k/4104), t+k);
            auto y6 = f(U[n] - y1*(8*k/27) + y2*(2*k) - y3*(3544*k/2565)
                        + y4*(1859*k/4104) - y5*(11*k/40), t+0.5*k);
            auto u1 = U[n] + y1*(25*k/216) + y3*(1408*k/2565) + y4*(2197*k/4104) - y5*(0.2*k);
            auto u2 = U[n] + y1*(16*k/135) + y3*(6656*k/12825) + y4*(28561*k/56430)
                      - y5*(9*k/50) + y6*(2*k/55);
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
std::unique_ptr<TimeIntegrator<dim>> createF45RK(int steps=100)
{
    return std::make_unique<F45embedRKsolver<dim>>(steps);
}


#endif