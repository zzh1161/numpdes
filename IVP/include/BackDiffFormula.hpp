#ifndef __BACK_DIFF_FORMULA_HPP__
#define __BACK_DIFF_FORMULA_HPP__

#include"TimeIntegrator.hpp"

template<int dim>
class BackDiffFormulaSolver : public TimeIntegrator<dim>
{
public:
    using Func = std::function<Vec<real,dim>(Vec<real,dim>,real)>;

    BackDiffFormulaSolver(int steps): steps_(steps){}
    BackDiffFormulaSolver(const BackDiffFormulaSolver &) = delete;
    BackDiffFormulaSolver& operator=(const BackDiffFormulaSolver &) = delete;
    ~BackDiffFormulaSolver() = default;

    Vec<real,dim> solve(const Vec<real,dim> &U0, Func f, real T, int p_=0) const
    {
        #define U  TimeIntegrator<dim>::result_
        std::vector<Vec<real,dim>>().swap(U);
        real k = 1.0*T/steps_;
        real epsilon = 0.000001;
        U.push_back(U0);
        for(int i=0; i<p_-1; i++){
            auto y1 = f(U[i], i*k);
            auto y2 = f(U[i]+y1*(0.5*k), i*k+0.5*k);
            auto y3 = f(U[i]+y2*(0.5*k), i*k+0.5*k);
            auto y4 = f(U[i]+y3*k, (i+1)*k);
            U.push_back(U[i] + (y1+y2*2+y3*2+y4)*(1.0*k/6));
        }
        if(p_ == 1){
            for(int n=1; n<=steps_; n++){
                auto u = U[n-1];
                auto g = [&](Vec<real,dim> &x){
                    return (f(x,n*k)*k + U[n-1]);
                };
                while(true){
                    u = g(u);
                    real err = norm(u-g(u), 0);
                    if(err < epsilon) break;
                }
                U.push_back(u);
            }
        }else if(p_ == 2){
            for(int n=2; n<=steps_; n++){
                auto u = U[n-1];
                auto g = [&](Vec<real,dim> &x){
                    return (f(x,n*k)*(2.0*k/3) + U[n-1]*(4.0/3) - U[n-2]*(1.0/3));
                };
                while(norm(u-g(u), 0)>=epsilon){
                    u = g(u);
                }
                U.push_back(u);
            }
        }else if(p_ == 3){
            for(int n=3; n<=steps_; n++){
                auto u = U[n-1];
                auto g = [&](Vec<real,dim> &x){
                    return (f(x,n*k)*(6.0*k/11) + U[n-1]*(18.0/11) - U[n-2]*(9.0/11)
                                                + U[n-3]*(2.0/11));
                };
                while(norm(u-g(u), 0)>=epsilon){
                    u = g(u);
                }
                U.push_back(u);
            }
        }else if(p_ == 4){
            for(int n=4; n<=steps_; n++){
                auto u = U[n-1];
                auto g = [&](Vec<real,dim> &x){
                    return (f(x,n*k)*(12.0*k/25) + U[n-1]*(48.0/25) - U[n-2]*(36.0/25)
                                                 + U[n-3]*(16.0/25) - U[n-4]*(3.0/25));
                };
                while(norm(u-g(u), 0)>=epsilon){
                    u = g(u);
                }
                U.push_back(u);
            }
        }else{
            throw std::runtime_error("No such p for BDF method!");
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
std::unique_ptr<TimeIntegrator<dim>> createBDF(int steps)
{
    return std::make_unique<BackDiffFormulaSolver<dim>>(steps=100);
}


#endif