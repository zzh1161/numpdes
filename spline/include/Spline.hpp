/**
 * @file Spline.hpp
 * @author zhang-zh
 * @brief A header file for Spline class.
 * @version 0.1
 * @date 2021-11-16
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef __SPLINE_HPP__
#define __SPLINE_HPP__

#include"Polynomial.hpp"
#include"InterpConditions.hpp"
#include<lapacke.h>

enum SplineType{ppForm, cardinalB};
enum BCType{complete, notAknot, periodic};

/**
 * @tparam Dim       Dimension of the spline.
 * @tparam Order     Order of the spline.
 * @tparam st        Type of the spline.
 * @tparam CoefType  Type of the polynomials' coefficients.
 * @tparam KnotType  Type of the knots.
 */
template<int Dim, int Order, SplineType st, class CoefType, class KnotType>
class Spline
{
public:
    std::vector<Polynomial<Order-1,CoefType>> spline;
    std::vector<KnotType> knots;
};

/**
 * A partial specialization to implement
 * ppForm spline.
 */
template<int Dim, int Order, class CoefType, class KnotType>
class Spline<Dim,Order,ppForm,CoefType,KnotType>
{
public:
    std::vector<Polynomial<Order-1,CoefType>> spline;
    std::vector<KnotType> knots;

public:
    /**
     * A method to generate an 
     * one-dimension ppForm spline.
     */
    template<int Ord>
    friend Spline<1,Ord,ppForm,double,double> interpolate(const InterpConditions &interp, BCType t);

    /**
     * Generate two-dimension ppForm splines.
     */
    template<int Ord>
    friend Spline<2,Ord,ppForm,Vec<double,2>,double> fitCurve(const std::vector<Vec<double,2>> &vec, BCType type);

    /**
     * Calculate at some x
     */
    template<class T>
    auto operator()(const T &x)
    {
        using Tx = decltype(spline[0](x));
        Tx res = static_cast<Tx>(0);
        for(int i=0; i<knots.size()-1; i++)
            if(knots[i]<=x && knots[i+1]>x){
                res = spline[i](x);
                break;
            }
        return res;
    }

    /**
     * Redirect the output of a ppForm spline.
     */
    friend std::ostream& operator<<(std::ostream &os, const Spline<Dim,Order,ppForm,CoefType,KnotType> &sp)
    {
        int n = sp.spline.size();
        os << "====== It's a ppForm spline. ======" << std::endl;
        os << "The knots are: " << std::endl;
        for(auto &index : sp.knots)
            os << index << " ";
        os << std::endl;
        os << "The spline is: " << std::endl;
        os << sp.spline[0] << ".*(x>=" << sp.knots[0]
        << " & x<" << sp.knots[1] << ") ..." << std::endl;
        for(int i=1; i<n-1; i++)
            os <<"+"<< sp.spline[i] << ".*(x>=" << sp.knots[i]
            << " & x<" << sp.knots[i+1] << ") ..." << std::endl;
        os << "+" << sp.spline[n-1] << ".*(x>=" << sp.knots[n-1]
        << " & x<" << sp.knots[n] << ");" << std::endl;
        os << "======================================" << std::endl;
        return os;
    }

    template<int dim, int ord>
    friend void output(const Spline<dim,ord,ppForm,Vec<double,dim>,double> &p);
};

template<int dim, int ord>
void output(const Spline<dim,ord,ppForm,Vec<double,dim>,double> &sp)
{
    int n = sp.spline.size();
#define os std::cout
    for(int j=0; j<dim; j++){
        os << "The spline of No." << j+1 << " is:" << std::endl;
        output(sp.spline[0],j);
        os << ".*(x>=" << sp.knots[0] << " & x<" << sp.knots[1] << ") ..." << std::endl;
        for(int i=1; i<n-1; i++){
            os << "+";
            output(sp.spline[i],j);
            os << ".*(x>=" << sp.knots[i] << " & x<" << sp.knots[i+1] << ") ..." << std::endl;
        }
        os << "+";
        output(sp.spline[n-1],j);
        os << ".*(x>=" << sp.knots[n-1] << " & x<" << sp.knots[n] << ");" << std::endl;
        os << std::endl;
    }
#undef os
}

/**
 * A partial specialization to implement
 * cardinalB spline.
 */
template<int Dim, int Order, class CoefType, class KnotType>
class Spline<Dim,Order,cardinalB,CoefType,KnotType>
{
public:
    std::vector<Polynomial<Order,CoefType>> spline;
    std::vector<KnotType> knots;
    int begindex; // Beginning of the knots
    int size;  // Number of knots;

public:

    /** 
     * addition and substraction for cardinal-B spline
     */
#define OPERATION_ADD_SUB(OpNm,Op)                             \
    auto OpNm(const Spline<Dim,Order,cardinalB,CoefType,KnotType> &p) const \
    {                                                          \
        int k = p.begindex - begindex;                         \
        assert(k<=size && k>=-p.size);                         \
        Spline<Dim,Order,cardinalB,CoefType,KnotType> res;     \
        res.begindex = (k>0?begindex:p.begindex);              \
        if(k>=0){                                              \
            if(p.size+k >= size){                              \
                res.size = p.size + k;                         \
                res.knots.resize(res.size);                    \
                res.spline.resize(res.size-1);                 \
                for(int i=0; i<size; res.knots[i]=knots[i],i++);\
                for(int i=size; i<p.size+k; res.knots[i]=p.knots[i-k],i++);\
                for(int i=0; i<k; i++)                         \
                    res.spline[i] = spline[i];                 \
                for(int i=k; i<size-1; i++)                    \
                    res.spline[i]=spline[i] Op p.spline[i-k];  \
                for(int i=size-1; i<p.size+k-1; i++)           \
                    res.spline[i]=Op p.spline[i-k];            \
            }                                                  \
            else{                                              \
                res.size = size;                               \
                res.knots.resize(res.size);                    \
                res.spline.resize(res.size-1);                 \
                for(int i=0; i<size; res.knots[i]=knots[i],i++);\
                for(int i=0; i<k; i++)                         \
                    res.spline[i]=spline[i];                   \
                for(int i=k; i<k+p.size-1; i++)                \
                    res.spline[i]=spline[i] Op p.spline[i-k];  \
                for(int i=k-1+p.size; i<size-1; i++)           \
                    res.spline[i]=spline[i];                   \
            }                                                  \
        }                                                      \
        else{                                                  \
            if(p.size+k <= size){                              \
                res.size = size-k;                             \
                res.knots.resize(res.size);                    \
                res.spline.resize(res.size-1);                 \
                for(int i=0; i<p.size; res.knots[i]=p.knots[i],i++);\
                for(int i=p.size; i<res.size; res.knots[i]=knots[i],i++);\
                for(int i=0; i<-k; i++)                        \
                    res.spline[i] =Op p.spline[i];             \
                for(int i=-k; i<p.size-1; i++)                 \
                    res.spline[i]=spline[i+k] Op p.spline[i];  \
                for(int i=p.size-1; i<size-k-1; i++)           \
                    res.spline[i] = spline[i+k];               \
            }                                                  \
            else{                                              \
                res.size=p.size;                               \
                res.knots.resize(res.size);                    \
                res.spline.resize(res.size-1);                 \
                for(int i=0; i<res.size; res.knots[i]=p.knots[i],i++);\
                for(int i=0; i<-k; i++)                        \
                    res.spline[i] = Op p.spline[i];            \
                for(int i=-k; i<size-k-1; i++)                 \
                    res.spline[i]= spline[i+k] Op p.spline[i]; \
                for(int i=size-k-1; i<p.size-1; i++)           \
                    res.spline[i] = Op p.spline[i];            \
            }                                                  \
        }                                                      \
        return res;                                            \
    }

    OPERATION_ADD_SUB(operator+, +)
    OPERATION_ADD_SUB(operator-, -)
#undef OPERATION_ADD_SUB

    /**
     * Multiplied by a polynomial
     */
    template<int num>
    auto operator*(const Polynomial<num,CoefType> &p) const
    {
        const int ord = Order+num;
        Spline<Dim,ord,cardinalB,CoefType,KnotType> res;
        res.begindex = begindex;
        res.size = size;
        res.knots.resize(size);
        res.spline.resize(size-1);
        for(int i=0; i<size; i++)
            res.knots[i] = knots[i];
        for(int i=0; i<size-1; i++)
            res.spline[i] = spline[i]*p;
        return res;
    }

    /**
     * Multiplied by a number
     */
    template<class T>
    auto operator*(const T &p) const
    {
        using Tx = decltype(spline[0][0]*p);
        Spline<Dim,Order,cardinalB,Tx,KnotType> res;
        res.begindex = begindex;
        res.size = size;
        res.knots.resize(size);
        res.spline.resize(size-1);
        for(int i=0; i<res.size; i++)
            res.knots[i] = knots[i];
        for(int i=0; i<size-1; i++)
            res.spline[i] = spline[i]*p;
        return res;
    }

    /**
     * Calculate at some x
     */
    template<class T>
    auto operator()(const T &x)
    {
        using Tx = decltype(spline[0](x));
        Tx res = static_cast<Tx>(0);
        for(int i=0; i<knots.size()-1; i++)
            if(knots[i]<=x && knots[i+1]>x){
                res = spline[i](x);
                break;
            }
        return res;
    }

    /**
     * Generate a cardinalB spline B_i^n with recursion
     */
    template<int n>
    friend Spline<1,n,cardinalB,double,double> Generate(int i);

    /**
     * A method to interpolate a series of points
     * to get an one-dimension cardinalB spline.
     */
    template<int Ord>
    friend Spline<1,Ord,cardinalB,double,double> interpolate(const InterpConditions &interp);

    /**
     * Redirect the output of a cardinalB spline.
     */
    friend std::ostream& operator<<(std::ostream &os, const Spline<Dim,Order,cardinalB,CoefType,KnotType> &sp)
    {
        int n = sp.spline.size();
        os << "====== It's a cardinalB spline. ======" << std::endl;
        os << "The knots are: " << std::endl;
        for(auto &index : sp.knots)
            os << index << " ";
        os << std::endl;
        os << "Polynomials are: " << std::endl;
        os << sp.spline[0] << ".*(x>=" << sp.knots[0]
        << " & x<" << sp.knots[1] << ") ..." << std::endl;
        for(int i=1; i<n-1; i++)
            os <<"+"<< sp.spline[i] << ".*(x>=" << sp.knots[i]
            << " & x<" << sp.knots[i+1] << ") ..." << std::endl;
        os << "+" << sp.spline[n-1] << ".*(x>=" << sp.knots[n-1]
        << " & x<" << sp.knots[n] << ");" << std::endl;
        os << "======================================" << std::endl;
        return os;
    }
};

template<int n>
Spline<1,n,cardinalB,double,double> Generate(int i)
{
    Polynomial<1,double> a{(-i+1)*1.0/n, 1.0/n};
    Polynomial<1,double> b{(i+n)*1.0/n,-1.0/n};
    return Generate<n-1>(i)*a + Generate<n-1>(i+1)*b;
}
template<>
Spline<1,1,cardinalB,double,double> Generate<1>(int i)
{
    Spline<1,1,cardinalB,double,double> res;
    res.begindex = i-1;
    res.size = 3;
    res.knots.resize(3);
    res.spline.resize(2);
    for(int j=0; j<res.size; j++)
        res.knots[j] = i-1+j;
    res.spline[0] = Polynomial<1,double>{(-i+1)*1.0, 1};
    res.spline[1] = Polynomial<1,double>{(i+1)*1.0, -1};
    return res;
}

////////////////////////////////////////////////////////
template<int Ord>
Spline<1,Ord,ppForm,double,double> interpolate(const InterpConditions &interp, BCType t);
/**
 * One-dimension interpolation when order=4
 */
template<>
Spline<1,4,ppForm,double,double> interpolate(const InterpConditions &interp, BCType t)
{
#define FV interp.FunctionValue
#define IP interp.InterpPoints
#define N interp.num
#define mu(i) (IP[i]-IP[i-1])/(IP[i+1]-IP[i-1])
#define lam(i) (IP[i+1]-IP[i])/(IP[i+1]-IP[i-1])

    assert(N >= 2);
    Spline<1,4,ppForm,double,double> Res;
    Res.knots.resize(N);
    Res.spline.resize(N-1);
    double B[N]; // right-hand sides
    for(int i=1; i<N-1; i++){
        B[i] = 6*(-(FV[i][0]-FV[i-1][0])/(IP[i]-IP[i-1])
                + (FV[i+1][0]-FV[i][0])/(IP[i+1]-IP[i]))
                /(IP[i+1]-IP[i-1]);
    }
    double A[N*N]; // Matrix of coefficient
    for(int i=1; i<N-1; i++){
        A[i*N+i-1] = mu(i);
        A[i*N+i] = 2;
        A[i*N+i+1] = lam(i);
        for(int j=0; j<N; j++){
            if(j!=i-1 && j!=i && j!=i+1)
                A[i*N+j] = 0;
        }
    }
    Polynomial<1,double> poly_1[N]; // x-x_i
    Polynomial<2,double> poly_2[N]; // (x-x_i)^2
    Polynomial<3,double> poly_3[N]; // (x-x_i)^3
    for(int i=0; i<N; i++){
        poly_1[i].set({-IP[i],1});
        poly_2[i].set({IP[i]*IP[i],-2*IP[i],1});
        poly_3[i].set({-pow(IP[i],3),3*pow(IP[i],2),-3*IP[i],1});
        Res.knots[i] = IP[i];
    }
#define GENERATE_SPLINE                                                \
    int ipiv[N];                                                       \
    LAPACKE_dgesv(LAPACK_ROW_MAJOR,N,1,A,N,ipiv,B,1);                  \
    for(int i=0; i<N-1; i++){                                          \
        double s = (FV[i+1][0]-FV[i][0])/(IP[i+1]-IP[i])               \
                    - (B[i+1]+2*B[i])*(IP[i+1]-IP[i])/6;               \
        Res.spline[i] = poly_3[i]*((B[i+1]-B[i])/(6*(IP[i+1]-IP[i])))  \
                        + poly_2[i]*(B[i]/2) +                         \
                        poly_1[i]*s + Polynomial<0,double>{FV[i][0]};  \
    }

    if(t == complete){
        assert(FV[0].size()>=2 && FV.back().size()>=2);
        B[0] = 6*((FV[1][0]-FV[0][0])/(IP[1]-IP[0]) 
                 - FV[0][1])/(IP[1]-IP[0]);
        B[N-1] = 6*(FV[N-1][1] - 
                    (FV[N-1][0]-FV[N-2][0])/
                    (IP[N-1]-IP[N-2]))/
                    (IP[N-1]-IP[N-2]);
        A[0] = 2; A[1] = 1;
        for(int j=2; j<N; j++)
            A[j] = 0;
        A[N*N-1] = 2; A[N*N-2] = 1;
        for(int j=0; j<N-2; j++)
            A[(N-1)*N+j] = 0;
        GENERATE_SPLINE
    }
    else if(t == notAknot){
        assert(N>=4);
        B[0] = 0; B[N-1] = 0;
        A[0]=IP[2]-IP[1]; A[1]=IP[0]-IP[2]; A[2]=IP[1]-IP[0];
        for(int j=3; j<N; j++)
            A[j]=0;
        A[N*N-1] = IP[N-2]-IP[N-3];
        A[N*N-2] = IP[N-3]-IP[N-1];
        A[N*N-3] = IP[N-1]-IP[N-2];
        for(int j=0; j<N-3; j++)
            A[N*(N-1)+j] = 0;
        GENERATE_SPLINE
    }
    else{
        assert(N>=4);
        B[0] = 0;
        B[N-1] = 6*(FV[1][0]-FV[0][0])/(IP[1]-IP[0]) - 
                 6*(FV[N-1][0]-FV[N-2][0])/(IP[N-1]-IP[N-2]);
        A[0] = 1; A[N-1] = -1;
        for(int i=1; i<N-1; A[i++]=0);
        A[N*(N-1)] = 2*(IP[1]-IP[0]);
        A[N*(N-1)+1] = IP[1]-IP[0];
        A[N*N-2] = IP[N-1]-IP[N-2];
        A[N*N-1] = 2*A[N*N-2];
        for(int j=2; j<N-2; A[N*N-N+j]=0,j++);
        GENERATE_SPLINE
    }

#undef GENERATE_SPLINE
#undef lam
#undef mu
#undef N
#undef IP
#undef FV
    return Res;
}
/**
 * One-dimension interpolation when order=2
 */
template<>
Spline<1,2,ppForm,double,double> interpolate(const InterpConditions &interp, BCType t)
{
#define FV interp.FunctionValue
#define IP interp.InterpPoints
#define N interp.num

    assert(N>=2);
    Spline<1,2,ppForm,double,double> Res;
    Res.knots.resize(N);
    Res.spline.resize(N-1);
    Polynomial<1,double> poly[N]; // x-x_i
    for(int i=0; i<N; i++){
        poly[i].set({-IP[i],1});
        Res.knots[i] = IP[i];
    }
    for(int j=0; j<N-1; j++)
        Res.spline[j] = poly[j]*((FV[j+1][0]-FV[j][0])/(IP[j+1]-IP[j]))
                        + FV[j][0];

#undef N
#undef IP
#undef FV
    return Res;
}

///////////////////////////////////////////////////////
template<int Ord>
Spline<2,Ord,ppForm,Vec<double,2>,double> fitCurve(const std::vector<Vec<double,2>> &vec, BCType type);
/**
 * Two-dimension fitCurve when order=4
 * 
 * Here we claim that if it's complete cubic spline,
 * the last two elements of vec are the fixed
 * derivatives at t_1 and t_N.
 */
template<>
Spline<2,4,ppForm,Vec<double,2>,double> fitCurve(const std::vector<Vec<double,2>> &vec, BCType type)
{
#define mu(i) (t[i]-t[i-1])/(t[i+1]-t[i-1])
#define lam(i) (t[i+1]-t[i])/(t[i+1]-t[i-1])
#define Ve Vec<double,2>
#define GENERATE_t                           \
    double t[N];                             \
    t[0]=0; Res.knots[0]=0;                  \
    for(int i=1; i<N; i++){                  \
        t[i] = t[i-1]+norm(vec[i]-vec[i-1]); \
        Res.knots[i] = t[i];                 \
    }
#define GENERATE_A                           \
    double A[N*N];                           \
    for(int i=1; i<N-1; i++){                \
        A[i*N+i-1] = mu(i);                  \
        A[i*N+i] = 2;                        \
        A[i*N+i+1] = lam(i);                 \
        for(int j=0; j<N; j++){              \
            if(j!=i-1 && j!=i && j!=i+1)     \
                A[i*N+j] = 0;                \
        }                                    \
    }
#define GENERATE_Spline                                \
    int ipiv[N];                                       \
    LAPACKE_dgesv(LAPACK_ROW_MAJOR,N,1,A,N,ipiv,b1,1); \
    LAPACKE_dgesv(LAPACK_ROW_MAJOR,N,1,A_,N,ipiv,b2,1);\
    Polynomial<1,double> poly_1[N];                    \
    Polynomial<2,double> poly_2[N];                    \
    Polynomial<3,double> poly_3[N];                    \
    for(int i=0; i<N; i++){                            \
        poly_1[i].set({-t[i],1});                      \
        poly_2[i].set({t[i]*t[i],-2*t[i],1});          \
        poly_3[i].set({-t[i]*t[i]*t[i],3*t[i]*t[i],-3*t[i],1}); \
    }                                                  \
    for(int i=0; i<N-1; i++){                          \
        Res.spline[i] = poly_3[i]*((Ve{b1[i+1],b2[i+1]}-Ve{b1[i],b2[i]})/((t[i+1]-t[i])*6))                                 \
                        + poly_2[i]*(Ve{b1[i],b2[i]}/2) +                                                                   \
                        poly_1[i]*((vec[i+1]-vec[i])/(t[i+1]-t[i])-(Ve{b1[i+1],b2[i+1]}+Ve{b1[i],b2[i]}*2)*(t[i+1]-t[i])/6) \
                        + Polynomial<0,Ve>{vec[i]};    \
    }

    Spline<2,4,ppForm,Vec<double,2>,double> Res;
    if(type == complete){
        int N = vec.size()-2;
        assert(N >= 2);
        Res.knots.resize(N);
        Res.spline.resize(N-1);
        GENERATE_t
        GENERATE_A
        A[0] = 2; A[1] = 1;
        for(int j=2; j<N; j++)
            A[j] = 0;
        A[N*N-1] = 2; A[N*N-2] = 1;
        for(int j=0; j<N-2; j++)
            A[(N-1)*N+j] = 0;
        double A_[N*N]; // As a copy of A
        for(int i=0; i<N*N; A_[i]=A[i],i++);
        Vec<double,2> B[N]; // right-hand sides
        for(int i=1; i<N-1; i++){
            B[i] = (((vec[i+1]-vec[i])/(t[i+1]-t[i])-(vec[i]-vec[i-1])/(t[i]-t[i-1]))/(t[i+1]-t[i-1]))*6;
        }
        B[0] = (((vec[1]-vec[0])/(t[1]-t[0])-vec[N])/(t[1]-t[0]))*6;
        B[N-1] = ((vec[N+1]-(vec[N-1]-vec[N-2])/(t[N-1]-t[N-2]))/(t[N-1]-t[N-2]))*6;
        double b1[N], b2[N]; // the 1st and 2nd part of B
        for(int i=0; i<N; i++){
            b1[i] = B[i][0];
            b2[i] = B[i][1];
        }
        GENERATE_Spline
    }
    else if(type == notAknot){
        int N = vec.size();
        assert(N >= 4);
        Res.knots.resize(N);
        Res.spline.resize(N-1);
        GENERATE_t
        GENERATE_A
        A[0]=t[2]-t[1]; A[1]=t[0]-t[2]; A[2]=t[1]-t[0];
        for(int j=3; j<N; j++)
            A[j]=0;
        A[N*N-1] = t[N-2]-t[N-3];
        A[N*N-2] = t[N-3]-t[N-1];
        A[N*N-3] = t[N-1]-t[N-2];
        for(int j=0; j<N-3; j++)
            A[N*(N-1)+j] = 0;
        double A_[N*N]; // As a copy of A
        for(int i=0; i<N*N; A_[i]=A[i],i++);
        Vec<double,2> B[N]; // right-hand sides
        for(int i=1; i<N-1; i++){
            B[i] = (((vec[i+1]-vec[i])/(t[i+1]-t[i])-(vec[i]-vec[i-1])/(t[i]-t[i-1]))/(t[i+1]-t[i-1]))*6;
        }
        B[0] = Ve{0,0}; B[N-1] = Ve{0,0};
        double b1[N], b2[N]; // the 1st and 2nd part of B
        for(int i=0; i<N; i++){
            b1[i] = B[i][0];
            b2[i] = B[i][1];
        }
        GENERATE_Spline
    }
    else{
        int N = vec.size();
        assert(N >= 4);
        Res.knots.resize(N);
        Res.spline.resize(N-1);
        GENERATE_t
        GENERATE_A
        A[0] = 1; A[N-1] = -1;
        for(int i=1; i<N-1; A[i++]=0);
        A[N*(N-1)] = 2*(t[1]-t[0]);
        A[N*(N-1)+1] = t[1]-t[0];
        A[N*N-2] = t[N-1]-t[N-2];
        A[N*N-1] = 2*A[N*N-2];
        for(int j=2; j<N-2; A[N*N-N+j]=0,j++);
        double A_[N*N]; // As a copy of A
        for(int i=0; i<N*N; A_[i]=A[i],i++);
        Vec<double,2> B[N]; // right-hand sides
        for(int i=1; i<N-1; i++){
            B[i] = (((vec[i+1]-vec[i])/(t[i+1]-t[i])-(vec[i]-vec[i-1])/(t[i]-t[i-1]))/(t[i+1]-t[i-1]))*6;
        }
        B[0] = Ve{0,0};
        B[N-1] = ((vec[1]-vec[0])/(t[1]-t[0]))*6- 
                 ((vec[N-1]-vec[N-2])/(t[N-1]-t[N-2]))*6;
        double b1[N], b2[N]; // the 1st and 2nd part of B
        for(int i=0; i<N; i++){
            b1[i] = B[i][0];
            b2[i] = B[i][1];
        }
        GENERATE_Spline
    }

#undef GENERATE_Spline
#undef GENERATE_A
#undef GENERATE_t
#undef Ve
#undef lam
#undef mu
    return Res;
}
/**
 * Two-dimension fitCurve when order=2
 */
template<>
Spline<2,2,ppForm,Vec<double,2>,double> fitCurve(const std::vector<Vec<double,2>> &vec, BCType type)
{
    int N = vec.size();
    assert(N >= 2);
    Spline<2,2,ppForm,Vec<double,2>,double> Res;
    Res.knots.resize(N);
    Res.spline.resize(N-1);
    double t[N];                             
    t[0]=0; Res.knots[0]=0;                  
    for(int i=1; i<N; i++){                  
        t[i] = t[i-1]+norm(vec[i]-vec[i-1]); 
        Res.knots[i] = t[i];                 
    }
    Polynomial<1,double> poly[N];
    for(int i=0; i<N; i++)
        poly[i].set({-t[i],1});
    for(int j=0; j<N-1; j++)
        Res.spline[j] = poly[j]*((vec[j+1]-vec[j])/(t[j+1]-t[j]))+vec[j];
    return Res;
}

///////////////////////////////////////////////////////
template<int Ord>
Spline<1,Ord,cardinalB,double,double> interpolate(const InterpConditions &interp);
/**
 * Generate cubic B-splines, i.e. Corollary 4.58
 */
template<>
Spline<1,3,cardinalB,double,double> interpolate(const InterpConditions &interp)
{
#define N interp.num
#define x interp.InterpPoints
#define y interp.FunctionValue
    // Here we claim that x[i] = t_i with x[i+1]-x[i]=1
    // and y[i][0] = f(x[i])
    assert(N>=2 && y[0].size()>1 && y[N-1].size()>1);
    for(int i=0; i<N-1; i++)
        assert(x[i+1]-x[i] == 1);
    double B[N];  // right-hand sides
    for(int i=1; i<N-1; i++)
        B[i] = 6*y[i][0];
    B[0] = 6*y[0][0] + 2*y[0][1];
    B[N-1] = 6*y[N-1][0] - 2*y[N-1][1];
    double A[N*N]={0}; // Coefficient matrix
    A[0]=4; A[1]=2; 
    A[N*N-1]=4; A[N*N-2]=2;
    for(int i=1; i<N-1; i++){
        A[i*N+i-1]=1; A[i*N+i]=4; A[i*N+i+1]=1;
    }
    int ipiv[N];
    LAPACKE_dgesv(LAPACK_ROW_MAJOR,N,1,A,N,ipiv,B,1);
    auto res = Generate<3>(x[0]-2)*(B[1]-2*y[0][1]);
    for(int i=0; i<N; i++)
        res = res+Generate<3>(x[i]-1)*B[i];
    res = res+Generate<3>(x.back())*(B[N-2]+2*y[N-1][1]);
    return res;
#undef y
#undef x
#undef N
}
/**
 * Generate quadratic B-splines, i.e. Corollary 4.59
 */
template<>
Spline<1,2,cardinalB,double,double> interpolate(const InterpConditions &interp)
{
#define N interp.num
#define x interp.InterpPoints
#define y interp.FunctionValue
    // Here we claim that y[i][0] = f(i+3/2)
    // and y[0][1] = f(t_1-1/2), y[N-1][1] = f(t_{N}-1/2)
    // and x[i] = t_i+1/2
    assert(N>=2 && y[0].size()>1 && y[N-1].size()>1);
    for(int i=0; i<N-1; i++)
        assert(x[i+1]-x[i] == 1);
    double B[N];  // right-hand sides
    for(int i=1; i<N-1; i++)
        B[i] = 8*y[i][0];
    B[0]=8*y[0][0]-2*y[0][1];
    B[N-1] = 8*y[N-1][0]-2*y[N-1][1];
    double A[N*N] = {0};
    A[0]=5; A[1]=1;
    A[N*N-1]=5; A[N*N-2]=1;
    for(int i=1; i<N-1; i++){
        A[i*N+i-1]=1; A[i*N+i]=6; A[i*N+i+1]=1;
    }
    int ipiv[N];
    LAPACKE_dgesv(LAPACK_ROW_MAJOR,N,1,A,N,ipiv,B,1);
    auto res = Generate<2>(x[0]-1)*(2*y[0][1]-B[0]);
    for(int i=0; i<N; i++)
        res = res+Generate<2>(x[i])*B[i];
    res = res+Generate<2>(x.back()+1)*(2*y[N-1][1]-B[N-1]);
    return res;

#undef y
#undef x
#undef N
}

#endif