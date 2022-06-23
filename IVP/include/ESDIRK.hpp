#ifndef __ESDIRK_HPP__
#define __ESDIRK_HPP__

#include"TimeIntegrator.hpp"

template<int dim>
class ESDIRKsolver : public TimeIntegrator<dim>
{
public:
    using Func = std::function<Vec<real,dim>(Vec<real,dim>,real)>;

    ESDIRKsolver(int steps): steps_(steps){}
    ESDIRKsolver(const ESDIRKsolver &) = delete;
    ESDIRKsolver& operator=(const ESDIRKsolver &) = delete;
    ~ESDIRKsolver() = default;

    Vec<real,dim> solve(const Vec<real,dim> &U0, Func f, real T, int p_=0) const
    {
        #define U TimeIntegrator<dim>::result_
        std::vector<Vec<real,dim>>().swap(U);
        real k = 1.0*T/steps_;
        U.push_back(U0);
        auto epsilon = 0.000000001;
        for(int n=0; n<steps_; n++){
            auto y1 = f(U[n], n*k);
            // calculate y2
            auto y2 = y1;
            while(true){
                auto temp = f(U[n]+y1*(0.25*k)+y2*(0.25*k), n*k+0.5*k);
                if(norm(y2-temp, 0) < epsilon) break;
                else y2 = temp;
            }
            // calculate y3
            auto y3 = y2;
            while(true){
                auto temp = f(U[n]+y1*(8611*k/62500)-y2*(1743*k/31250)+y3*(0.25*k), n*k+(83*k/250));
                if(norm(y3-temp, 0) < epsilon) break;
                else y3 = temp;
            }
            // calculate y4
            auto y4 = y3;
            while(true){
                auto temp = f(U[n]+y1*(5012029*k/34652500)-y2*(654441*k/2922500)
                              +y3*(174375*k/388108)+y4*(0.25*k), n*k+(31*k/50));
                if(norm(y4-temp, 0) < epsilon) break;
                else y4 = temp;
            }
            // calculate y5
            auto y5 = y4;
            while(true){
                auto temp = f(U[n]+y1*(15267082809*k/155376265600)-y2*(71443401*k/120774400)
                              +y3*(730878875*k/902184768)+y4*(2285395*k/8070912)+y5*(0.25*k),
                              n*k+(17*k/20));
                if(norm(y5-temp, 0) < epsilon) break;
                else y5 = temp;
            }
            // calculate y6
            auto y6 = y5;
            while(true){
                auto temp = f(U[n]+y1*(82889*k/524892)+y3*(15625*k/83664)+y4*(69875*k/102672)
                              -y5*(2260*k/8211)+y6*(0.25*k), n*k+k);
                if(norm(y6-temp, 0) < epsilon) break;
                else y6 = temp;
            }
            U.push_back(U[n]+y1*(82889*k/524892)+y3*(15625*k/83664)+y4*(69875*k/102672)
                        -y5*(2260*k/8211)+y6*(0.25*k));
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
std::unique_ptr<TimeIntegrator<dim>> createESDIRK(int steps=100)
{
    return std::make_unique<ESDIRKsolver<dim>>(steps);
}


#endif