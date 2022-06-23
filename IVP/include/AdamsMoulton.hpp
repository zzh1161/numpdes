#ifndef __ADAMSMOULTON_HPP__
#define __ADAMSMOULTON_HPP__

#include"TimeIntegrator.hpp"

template<int dim>
class AdamsMoultonSolver : public TimeIntegrator<dim>
{
public:
    using Func = std::function<Vec<real,dim>(Vec<real,dim>,real)>;

    AdamsMoultonSolver(int steps): steps_(steps){}
    AdamsMoultonSolver(const AdamsMoultonSolver &) = delete;
    AdamsMoultonSolver& operator=(const AdamsMoultonSolver &) = delete;
    ~AdamsMoultonSolver() = default;

    Vec<real,dim> solve(const Vec<real,dim> &U0, Func f, real T, int p_=0) const
    {
        std::vector<Vec<real,dim>>().swap(this->result_);
        real k = 1.0*T/steps_;
        // For Adams-Moulton method, actually the s we need is p-1,
        // but we need keep at least p values, because we need 
        // Adams-Bashforth method to carry out predictor-corrector process.
        Vec<real,dim> init_U[p_]; init_U[0] = U0; this->result_.push_back(U0);
        for(int i=0; i<p_-1; i++){
            auto y1 = f(init_U[i], i*k);
            auto y2 = f(init_U[i]+y1*(0.5*k), i*k+0.5*k);
            auto y3 = f(init_U[i]+y2*(0.5*k), i*k+0.5*k);
            auto y4 = f(init_U[i]+y3*k, (i+1)*k);
            init_U[i+1] = init_U[i] + (y1+y2*2+y3*2+y4)*(1.0*k/6);
            this->result_.push_back(init_U[i+1]);
        }
        std::deque<Vec<real,dim>> U(init_U, init_U+p_);
        if(p_ == 2){
            for(int n=1; n<steps_; n++){
                auto P = U.back() - f(U[0],(n-1)*k)*(0.5*k) + f(U[1],n*k)*(1.5*k);
                auto u = U.back() + (f(P,(n+1)*k)+f(U[1],n*k))*(0.5*k);
                this->result_.push_back(u);
                U.push_back(u);
                U.pop_front();
            }
        }else if(p_ == 3){
            for(int n=2; n<steps_; n++){
                auto P = U.back() + f(U[0],(n-2)*k)*(5.0*k/12) - f(U[1],(n-1)*k)*(4.0*k/3)
                                  + f(U[2],n*k)*(23.0*k/12);
                auto u = U.back() + f(P,(n+1)*k)*(5.0*k/12) + f(U[2],n*k)*(8.0*k/12)
                                  - f(U[1],(n-1)*k)*(1.0*k/12);
                this->result_.push_back(u);
                U.push_back(u);
                U.pop_front();
            }
        }else if(p_ == 4){
            for(int n=3; n<steps_; n++){
                auto P = U.back() - f(U[0],(n-3)*k)*(3.0*k/8) + f(U[1],(n-2)*k)*(37.0*k/24)
                                  - f(U[2],(n-1)*k)*(59.0*k/24) + f(U[3],n*k)*(55.0*k/24);
                auto u = U.back() + f(P,(n+1)*k)*(9.0*k/24) + f(U[3],n*k)*(19.0*k/24)
                                  - f(U[2],(n-1)*k)*(5.0*k/24) + f(U[1],(n-2)*k)*(1.0*k/24);
                this->result_.push_back(u);
                U.push_back(u);
                U.pop_front();
            }
        }else if(p_ == 5){
            for(int n=4; n<steps_; n++){
                auto u = U.back();
                auto g = [U,k,n,f](Vec<real,dim> &x){
                    return U.back() + f(x,(n+1)*k)*(251.0*k/720) + f(U[4],n*k)*(646.0*k/720)
                                    - f(U[3],(n-1)*k)*(264.0*k/720) + f(U[2],(n-2)*k)*(106.0*k/720)
                                    - f(U[1],(n-3)*k)*(19.0*k/720);
                };
                while(norm((u-g(u)),0)>=0.0000001){
                    u = g(u);
                }
                this->result_.push_back(u);
                U.push_back(u);
                U.pop_front();
            }
        }else{
            throw std::runtime_error("No such p for Adams-Moulton method!");
        }
        return U.back();
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
std::unique_ptr<TimeIntegrator<dim>> createAdamsMoulton(int steps=100)
{
    return std::make_unique<AdamsMoultonSolver<dim>>(steps);
}

#endif