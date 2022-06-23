#ifndef __ADAMSBASHFORTH_HPP__
#define __ADAMSBASHFORTH_HPP__

#include"TimeIntegrator.hpp"

template<int dim>
class AdamsBashforthSolver : public TimeIntegrator<dim>
{
public:
    using Func = std::function<Vec<real,dim>(Vec<real,dim>,real)>;

    AdamsBashforthSolver(int steps): steps_(steps){}
    AdamsBashforthSolver(const AdamsBashforthSolver &) = delete;
    AdamsBashforthSolver& operator=(const AdamsBashforthSolver &) = delete;
    ~AdamsBashforthSolver() = default;

    Vec<real,dim> solve(const Vec<real,dim> &U0, Func f, real T, int p_=0) const
    {
        std::vector<Vec<real,dim>>().swap(this->result_);
        real k = 1.0*T/steps_;
        if(p_ == 1){
            auto U = U0; this->result_.push_back(U);
            for(int n=0; n<steps_; n++){
                U = U + f(U,k*n)*k;
                this->result_.push_back(U);
            }
            return U;
        }else if(p_ == 2){
            // use classical Runge-Kutta method to get the initial values
            auto y1 = f(U0, 0);
            auto y2 = f(U0+y1*(0.5*k), 0.5*k);
            auto y3 = f(U0+y2*(0.5*k), 0.5*k);
            auto y4 = f(U0+y3*k, k);
            auto U1 = U0 + (y1+y2*2+y3*2+y4)*(1.0*k/6);
            std::deque<Vec<real,dim>> U{U0,U1};
            this->result_.push_back(U0); this->result_.push_back(U1);
            for(int n=1; n<steps_; n++){
                auto temp = U.back() - f(U.front(),(n-1)*k)*(0.5*k) + f(U.back(),n*k)*(1.5*k);
                this->result_.push_back(temp);
                U.push_back(temp);
                U.pop_front();
            }
            return U.back();
        }else if(p_ == 3){
            // use classical Runge-Kutta method to get the initial values
            Vec<real,dim> init_U[3]; init_U[0] = U0; this->result_.push_back(init_U[0]);
            for(int i=0; i<2; i++){
                auto y1 = f(init_U[i], i*k);
                auto y2 = f(init_U[i]+y1*(0.5*k), i*k+0.5*k);
                auto y3 = f(init_U[i]+y2*(0.5*k), i*k+0.5*k);
                auto y4 = f(init_U[i]+y3*k, (i+1)*k);
                init_U[i+1] = init_U[i] + (y1+y2*2+y3*2+y4)*(1.0*k/6);
                this->result_.push_back(init_U[i+1]);
            }
            std::deque<Vec<real,dim>> U(init_U, init_U+3);
            for(int n=2; n<steps_; n++){
                auto temp = U.back() + f(U[0],(n-2)*k)*(5.0*k/12) - f(U[1],(n-1)*k)*(4.0*k/3)
                                     + f(U[2],n*k)*(23.0*k/12);
                this->result_.push_back(temp);
                U.push_back(temp);
                U.pop_front();
            }
            return U.back();
        }else if(p_ == 4){
            // use classical Runge-Kutta method to get the initial values
            Vec<real,dim> init_U[4]; init_U[0] = U0; this->result_.push_back(init_U[0]);
            for(int i=0; i<3; i++){
                auto y1 = f(init_U[i], i*k);
                auto y2 = f(init_U[i]+y1*(0.5*k), i*k+0.5*k);
                auto y3 = f(init_U[i]+y2*(0.5*k), i*k+0.5*k);
                auto y4 = f(init_U[i]+y3*k, (i+1)*k);
                init_U[i+1] = init_U[i] + (y1+y2*2+y3*2+y4)*(1.0*k/6);
                this->result_.push_back(init_U[i+1]);
            }
            std::deque<Vec<real,dim>> U(init_U, init_U+4);
            for(int n=3; n<steps_; n++){
                auto temp = U.back() - f(U[0],(n-3)*k)*(3.0*k/8) + f(U[1],(n-2)*k)*(37.0*k/24)
                                     - f(U[2],(n-1)*k)*(59.0*k/24) + f(U[3],n*k)*(55.0*k/24);
                this->result_.push_back(temp);
                U.push_back(temp);
                U.pop_front();
            }
            return U.back();
        }else{
            throw std::runtime_error("No such p for Adams-Bashforth Method!");
        }
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
std::unique_ptr<TimeIntegrator<dim>> createAdamsBashf(int steps=100)
{
    return std::make_unique<AdamsBashforthSolver<dim>>(steps);
} 

#endif