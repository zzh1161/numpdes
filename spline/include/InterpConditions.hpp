/**
 * @file InterpConditions.hpp
 * @author zhang-zh
 * @brief A header file for Interpolation Conditions.
 * @version 0.1
 * @date 2021-11-15
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef __INTERPOLATION_CONDITIONS__
#define __INTERPOLATION_CONDITIONS__

#include<functional>
#include<iostream>
#include<vector>
#include<numeric>

struct InterpConditions
{
public:
    // Number of interpolation points
    int num;
    // Interpolation points
    std::vector<double> InterpPoints;
    // Derivated times at each point
    std::vector<int> DiffTimes;
    // Corresponding values
    std::vector< std::vector<double> > FunctionValue;

    InterpConditions(){};

    InterpConditions(int _n, double *a, int *b, std::function<double(double)> *funcs)
    {
        num = _n;
        InterpPoints.resize(num);
        DiffTimes.resize(num);
        for(int i=0; i<num; i++)
        {
            InterpPoints[i] = a[i];
            DiffTimes[i] = b[i];
        }
        FunctionValue.resize(num);
        for(int i=0; i<num; i++)
        {
            FunctionValue[i].resize(b[i]);
            for(int j=0; j<b[i]; j++)
            {
                FunctionValue[i][j] = funcs[j](a[i]);
            }
        }
    }

    InterpConditions(int _n, double *a, int *b, double *c)
    {
        num = _n;
        InterpPoints.resize(num);
        DiffTimes.resize(num);
        for(int i=0; i<num; i++)
        {
            InterpPoints[i] = a[i];
            DiffTimes[i] = b[i];
        }
        FunctionValue.resize(num);
        for(int i=0; i<num; i++)
        {
            FunctionValue[i].resize(b[i]);
        }
        int i=0, j=0;
        while(i<num)
        {
            for(int k=0; k<FunctionValue[i].size(); k++)
            {
                FunctionValue[i][k] = c[j];
                j++;
            }
            i++;
        }
    }

    int get_order()
    {
        return accumulate(DiffTimes.begin(), DiffTimes.end(), -1);
    }

    int get_which_x(int j)
    {
        for(int i=0; i<num; i++)
        {
            if(j<=DiffTimes[i])
                return i;
            else
                j -= DiffTimes[i];
        }
    }

    friend std::ostream& operator<<(std::ostream &os, const InterpConditions &I)
    {
        os << "Interpolation points are:" << std::endl;
        for(int i=0; i<I.num; i++)
            for(int j=0; j<I.DiffTimes[i]; j++)
                os << I.InterpPoints[i] << " ";
        os << std::endl;
        os << "Corresponding values are:" << std::endl;
        for(int i=0; i<I.num; i++)
            for(int j=0; j<I.FunctionValue[i].size(); j++)
                os << I.FunctionValue[i][j] << " ";
        os << std::endl;
        return os;
    }    
};



#endif