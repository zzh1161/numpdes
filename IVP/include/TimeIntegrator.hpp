#ifndef __TIME_INTEGRATOR_HPP__
#define __TIME_INTEGRATOR_HPP__

#include<cassert>
#include<cmath>
#include<deque>
#include<iostream>
#include<vector>
#include<functional>
#include<memory>
#include"Info.hpp"
#include"Vec.h"

template<int dim>
class TimeIntegrator
{
public:
    using Func = std::function<Vec<real,dim>(Vec<real,dim>,real)>;

    TimeIntegrator(){}
    virtual ~TimeIntegrator() = default;

    void display() const
    {
        for(auto i=result_.begin(); i<result_.end(); i++){
            std::cout << static_cast<int>(i-result_.begin()) << ": "
                      << *i << std::endl;
        }
    }

    void disp_for_matplot() const
    {
        for(int i=0; i<dim; i++){
            std::cout << "u[" << i+1 << "] = [";
            for(auto &vec : result_){
                std::cout << vec[i] << ", ";
            }
            std::cout << "]" << std::endl;
        }
    }

    friend std::ostream& operator<<(std::ostream &os, TimeIntegrator<dim> &TI)
    {
        for(int i=0; i<dim; i++){
            os << "u[" << i+1 << "] = [";
            int anotherline = 0;
            for(auto &vec : TI.result_){
                os << vec[i] << ", ";
                anotherline++;
                if(anotherline >= 20){
                    os << "..." << std::endl;
                    anotherline = 0;
                }
            }
            os << "]" << std::endl;
        }
        return os;
    }

    virtual Vec<real,dim> solve(const Vec<real,dim> &U0, Func f, real T, int p_=0) const = 0;
    virtual real solution_error(const Vec<real,dim> &U0, Func f, real T, int p_=0) = 0;
    
protected:
    mutable std::vector<Vec<real,dim>> result_;
};



#endif