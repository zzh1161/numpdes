/**
 * @file Polynomial.hpp
 * @author zhang-zh
 * @brief A header file for class Polynomial with template.
 * @version 0.1
 * @date 2021-11-14
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef __POLYNOMIAL_HPP__
#define __POLYNOMIAL_HPP__

#include<iostream>
#include<algorithm>
#include<initializer_list>
#include<memory>
#include<vector>
#include<cassert>
#include<cmath>
#include"Vec.h"

template<int Order, class CoefType>
class Polynomial
{
protected:
    // Coefficients of the polynomial
    CoefType _coefs[Order+1];

public:
    /**
     * Constructors
     */

    Polynomial(){};

    Polynomial(std::initializer_list<CoefType> lst)
    {
        assert(lst.size()>=Order+1);
        auto j = lst.begin();
        for(int i=0; i<=Order; i++)
            _coefs[i] = *j++;
    }

    void set(std::initializer_list<CoefType> lst)
    {
        assert(lst.size()>=Order+1);
        auto j = lst.begin();
        for(int i=0; i<=Order; i++)
            _coefs[i] = *j++;
    }

    template<class T>
    explicit Polynomial(const Polynomial<Order,T> &p)
    {
        for(int i=0; i<=Order; i++)
            _coefs[i] = static_cast<CoefType>(p[i]);
    }

    template<class T>
    void operator=(const Polynomial<Order,T> &p)
    {
        for(int i=0; i<=Order; i++)
            _coefs[i] = static_cast<CoefType>(p[i]);
    }

    // // Move constructor
    // Polynomial(Polynomial<Order,CoefType> &&p)
    // {
    //     for(int i=0; i<=Order; i++)
    //         _coefs[i] = p[i];
    //     p._coefs = nullptr;
    //     std::cout << "Move Constructor." << std::endl;
    // }

    /**
     * Accessors
     */
    
    CoefType &operator[](int _d) { return _coefs[_d]; }

    const CoefType &operator[](int _d) const { return _coefs[_d]; }

    const CoefType *coeffis() const { return &_coefs[0]; }

    /**
     * Operations
     */

    // Evaluation at some x
    template<class T2>
    auto operator()(T2 _x)
    {
        using Tx = decltype(_coefs[0]*_x);
        Tx res;
        if(Order>=1){
            res = static_cast<Tx>(_coefs[0]);
            for(int i=1; i<=Order; i++)
                res = res + _coefs[i]*pow(_x, i);
        }
        else 
            res = static_cast<Tx>(_coefs[0]);
        return res;
    }

    // Addition and Subtraction
#define OPERATION_ADD_SUB_OP(OpNm, Op)               \
    template<int num, class T>                       \
    auto OpNm (const Polynomial<num,T> &p) const {   \
        using Tx = decltype(_coefs[0] Op p[0]);      \
        const int ord = (Order>num?Order:num);       \
        Polynomial<ord,Tx> res;                      \
        for(int i=0; i<=ord; i++){                   \
            if(i<=Order && i<=num)                   \
                res[i] = _coefs[i] Op p[i];          \
            else if(i>Order)                         \
                res[i] = static_cast<Tx>(p[i]);      \
            else                                     \
                res[i] = static_cast<Tx>(_coefs[i]); \
        }                                            \
        return res;                                  \
    }

    OPERATION_ADD_SUB_OP(operator+, +)
    OPERATION_ADD_SUB_OP(operator-, -)   
#undef OPERATION_ADD_SUB_OP

    // Friend operators
#define FRIEND_ADD_SUB_OP(OpNm,Op)                         \
    friend auto OpNm(const Polynomial<Order,CoefType> &p1, \
                     const CoefType &p2) {                 \
        Polynomial<Order,CoefType> res;                    \
        res[0] = p1[0] Op p2;                              \
        for(int i=1; i<=Order; i++)                        \
            res[i] = p1[i];                                \
        return res;                                        \
    }

    FRIEND_ADD_SUB_OP(operator+, +)
    FRIEND_ADD_SUB_OP(operator-, -)
#undef FRIEND_ADD_SUB_OP

    // Multiplication by a polynomial
    template<int num, class T>
    auto operator*(const Polynomial<num,T> &p) const 
    {
        using Tx = decltype(_coefs[0]*p[0]);
        const int ord = Order+num;
        Polynomial<ord,Tx> res;
        res[0] = _coefs[0]*p[0];
        // for(int i=1; i<=ord; res[i++]=0);
        for(int i=1; i<=ord; i++){
            int k;
            for(k=0; k<=Order; k++){
                if(i-k>=0 && i-k<=num){
                    res[i] = _coefs[k]*p[i-k];
                    break;
                }
            }
            for(int j=k+1; j<=Order; j++){
                if(i-j>=0 && i-j<=num)
                    res[i] = res[i] + _coefs[j]*p[i-j];
            }
        }
        return res;
    }

    // Mutiplication by a number
    template<class T>
    friend auto operator*(const Polynomial<Order,CoefType> &p1,
                          const T &p2)
    {
        using Tx = decltype(p2*p1[0]);
        Polynomial<Order,Tx> res;
        for(int i=0; i<=Order; i++)
            res[i] = p2*p1[i];
        return res;
    }

    // Derivation
    auto diff() const
    {
        assert(Order>0);
        const int ord = Order-1;
        using Tx = decltype(_coefs[0]*2);
        Polynomial<ord,Tx> res;
        for(int i=0; i<=ord; i++){
            res[i] = _coefs[i+1]*(i+1);
        }
        return res;
    }

    friend std::ostream& operator<<(std::ostream &os, const Polynomial<Order,CoefType> &p)
    {
        os << "( " << p[0];
        for(int i=1; i<=Order; i++)
            os << " + " << p[i] << "*x.^" << i;
        os << " )";
        return os;
    }

    template<int ord, int n>
    friend void output(const Polynomial<ord,Vec<double,n>> &p);
};

template<int ord, int n>
void output(const Polynomial<ord,Vec<double,n>> &p,int i)
{
    assert(i<n);
#define os std::cout
    os << "( " << p[0][i];
    for(int j=1; j<=ord; j++)
        os << " + " << p[j][i] << "*x.^" << j;
    os << " )";
#undef os
}

#define POSI_NEG(OpNm,Op)                                     \
    template<int ord, class T>                                \
    inline Polynomial<ord,T> OpNm(const Polynomial<ord,T> &p) \
    {                                                         \
        Polynomial<ord,T> res;                                \
        for(int i=0; i<=ord; i++)                             \
            res[i] = Op p[i];                                 \
        return res;                                           \
    }

    POSI_NEG(operator+,+)
    POSI_NEG(operator-,-)
#undef POSI_NEG



#endif