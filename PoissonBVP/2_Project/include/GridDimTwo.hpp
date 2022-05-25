/**************************************************************
 > @Description    : Multigrid methods in two-dimensional grid
 > @Version        : 1.0
 > @Author         : zhang-zh
 > @Date           : 2022-04-22 15:15
 > @LastEditTime   : 2022-04-27 15:18
**************************************************************/
#ifndef __PRO2_GRID_DIMTWO_HPP__
#define __PRO2_GRID_DIMTWO_HPP__

#include"Grid.hpp"

using vec2D = std::vector<std::vector<real>>;

template<>
class Grid<2,Domin::regular>
{
public:
    explicit Grid(const int n, BoundaryCondition<2> &B, 
                        std::function<real(real,real)> f){
        assert((n > 0) && (n&(n-1))==0 && n>=MINGrid);
        real h         = 1.0/n;
        size_          = n;
        BoundCondi_    = B;
        RightHandFunc_ = f;
        points_.resize(n+1);
        RightHandSide_.resize(n+1);
        for(auto &i : points_)        i.resize(n+1);
        for(auto &j : RightHandSide_) j.resize(n+1);
        for(int i=1; i<n; i++)
            for(int j=1; j<n; j++)
                // (-1 -1 4 -1 -1)/h^2
                RightHandSide_[i][j] = f(j*h, i*h);
        // (3a1/2h+b1)v[j][0] + (-2a1/h)v[j][1] + (a1/2h)v[j][2] = sigma
        RightHandSide_[0][0] = BoundCondi_.leftBound(0,0);
        RightHandSide_[0][n] = BoundCondi_.downBound(1,0);
        RightHandSide_[n][0] = BoundCondi_.upperBound(0,1);
        RightHandSide_[n][n] = BoundCondi_.rightBound(1,1);
        for(int i=1; i<n; i++){
            RightHandSide_[0][i] = BoundCondi_.downBound(i*h,0);
            RightHandSide_[n][i] = BoundCondi_.upperBound(i*h,1);
            RightHandSide_[i][0] = BoundCondi_.leftBound(0,i*h);
            RightHandSide_[i][n] = BoundCondi_.rightBound(1,i*h);
        }
    }

private:
    /**
     * @param  w       weighted Jacobi's coefficient
     * @param  Rest    restriction operator
     * @param  Inter   interpolation operator
     */
    vec2D V_cycle(real w, vec2D v, vec2D f, int nu1, int nu2, 
                Operators::Restriction Rest, Operators::Interpolation Inter);
    vec2D FMG(real w, vec2D f, int nu1, int nu2, 
                Operators::Restriction Rest, Operators::Interpolation Inter);

public:
    void Poisson_BVP_Multigrid(MultigridMethod method, real w = 1, int nu1 = 3, int nu2 = 3, 
                               Operators::Restriction Rest = Operators::fullWeighting,
                               Operators::Interpolation Inter = Operators::linear,
                               int max_iter = 100, real epsilon = 0.0001);

    real analytic_error(std::function<real(real,real)> f) const{
        vec2D err(size_+1); for(auto &i:err) i.resize(size_+1);
        for(int i=0; i<=size_; i++)
            for(int j=0; j<=size_; j++)
                err[i][j] = abs(points_[i][j] - f(j*1.0/size_,i*1.0/size_));
        auto maxl = std::max_element(err.begin(), err.end(),
                        [](const std::vector<real> &lh, const std::vector<real> &rh){
                            return *std::max_element(lh.begin(),lh.end()) <              
                                   *std::max_element(rh.begin(),rh.end());
                        });
        return *std::max_element(maxl->begin(), maxl->end());
    }

    real discrete_error() const{
        auto al = BoundCondi_.leftMixed.first;  auto bl = BoundCondi_.leftMixed.second;
        auto ar = BoundCondi_.rightMixed.first; auto br = BoundCondi_.rightMixed.second;
        auto au = BoundCondi_.upperMixed.first; auto bu = BoundCondi_.upperMixed.second;
        auto ad = BoundCondi_.downMixed.first;  auto bd = BoundCondi_.downMixed.second;
        int m = (size_+1)*(size_+1), n = size_+1;
        real h = 1.0/size_;
        std::vector<real> A(m*m);
        std::vector<real> v(m);
        for(int i=1; i<size_; i++)
            for(int j=1; j<size_; j++){
                int index  = (i*n+j)*m + i*n+j;
                A[index]   =  4*size_*size_;
                A[index-n] = -1*size_*size_;
                A[index+n] = -1*size_*size_;
                A[index+1] = -1*size_*size_;
                A[index-1] = -1*size_*size_;
                v[i*n+j] = RightHandFunc_(j*h,i*h);
            }
        for(int i=0; i<size_; i++){
            int index  = i*n*m + i*n;
            A[index]   = 1.5*al*size_+bl;
            A[index+1] = -2*al*size_;
            A[index+2] = 0.5*al*size_;
            v[i*n] = BoundCondi_.leftBound(0,i*h);
            index = ((n-1)*n+i)*(m+1);
            A[index]     = 1.5*au*size_+bu;
            A[index-n]   = -2*au*size_;
            A[index-2*n] = 0.5*au*size_;
            v[(n-1)*n+i] = BoundCondi_.upperBound(i*h,1);
        }
        for(int i=1; i<=size_; i++){
            int index = i*(m+1);
            A[index]     = 1.5*ad*size_+bd;
            A[index+n]   = -2*ad*size_;
            A[index+2*n] = 0.5*ad*size_;
            v[i] = BoundCondi_.downBound(i*h,0);
            index = (i*n+n-1)*(m+1);
            A[index]   = 1.5*ar*size_+br;
            A[index-1] = -2*ar*size_;
            A[index-2] = 0.5*ar*size_;
            v[i*n+n-1] = BoundCondi_.rightBound(1,i*h);
        }
        int ipiv[m];
        LAPACKE_dgesv(LAPACK_ROW_MAJOR,m,1,&A[0],m,ipiv,&v[0],1);
        for(int i=0; i<=size_; i++)
            for(int j=0; j<=size_; j++)
                v[i*n+j] = abs(v[i*n+j]-points_[i][j]);
        return *std::max_element(v.begin(), v.end());
    }

    void display() const{
        real h = 1.0/size_;
        std::cout << "x = [";
        for(int i=0; i<=size_; i++)
            for(int j=0; j<=size_; j++)
                std::cout << j*h << ", ";
        std::cout << "];" << std::endl;
        std::cout << "y = [";
        for(int i=0; i<=size_; i++)
            for(int j=0; j<=size_; j++)
                std::cout << i*h << ", ";
        std::cout << "];" << std::endl;
        std::cout << "z = [";
        for(int i=0; i<=size_; i++)
            for(int j=0; j<=size_; j++)
                std::cout << points_[i][j] << ", ";
        std::cout << "];" << std::endl;
    }

    void error_convergece(std::function<real(real,real)> f);

private:
    int size_;
    std::vector<std::vector<real>> points_;
    BoundaryCondition<2>           BoundCondi_;
    std::function<real(real,real)> RightHandFunc_;
    std::vector<std::vector<real>> RightHandSide_;
};

/////////////////////////////////////////////////////////////////////////////////////
/******************************* V-cycle & FMG *************************************/
/////////////////////////////////////////////////////////////////////////////////////

vec2D Grid<2>::V_cycle(real w, vec2D v, vec2D f, int nu1, int nu2, 
                Operators::Restriction Rest, Operators::Interpolation Inter){
    int n = v.size()-1;
    assert((n > 0) && (n&(n-1))==0 && n>=MINGrid);
    assert(v.size()==n+1 && f.size()==n+1);
    real h = 1.0/n;
    auto al = BoundCondi_.leftMixed.first;  auto bl = BoundCondi_.leftMixed.second;
    auto ar = BoundCondi_.rightMixed.first; auto br = BoundCondi_.rightMixed.second;
    auto au = BoundCondi_.upperMixed.first; auto bu = BoundCondi_.upperMixed.second;
    auto ad = BoundCondi_.downMixed.first;  auto bd = BoundCondi_.downMixed.second;
    // Relax nu1 times
    vec2D temp(n+1); for(auto &i:temp) i.resize(n+1);
    for(int times=0; times<nu1; times++){
        for(int i=0; i<n; i++){
            temp[i][0] = v[i][0]*(1-w) + v[i][1]*4*al*n*w/(3*al*n+2*bl) - 
                         v[i][2]*al*n*w/(3*al*n+2*bl) + f[i][0]*2*w/(3*al*n+2*bl);
            temp[n][i] = v[n][i]*(1-w) + v[n-1][i]*4*au*n*w/(3*au*n+2*bu) - 
                         v[n-2][i]*au*n*w/(3*au*n+2*bu) + f[n][i]*2*w/(3*au*n+2*bu);
        }
        for(int i=1; i<=n; i++){
            temp[0][i] = v[0][i]*(1-w) + v[1][i]*4*ad*n*w/(3*ad*n+2*bd) -
                         v[2][i]*ad*n*w/(3*ad*n+2*bd) + f[0][i]*2*w/(3*ad*n+2*bd);
            temp[i][n] = v[i][n]*(1-w) + v[i][n-1]*4*ar*n*w/(3*ar*n+2*br) -
                         v[i][n-2]*ar*n*w/(3*ar*n+2*br) + f[i][n]*2*w/(3*ar*n+2*br);
        }
        for(int i=1; i<n; i++)
            for(int j=1; j<n; j++){
                temp[i][j] = v[i][j]*(1-w) + 0.25*w*v[i-1][j] + 0.25*w*v[i+1][j] +
                             0.25*w*v[i][j-1] + 0.25*w*v[i][j+1] + 0.25*h*h*w*f[i][j];
            }
        v = temp;
    }
    // Recursion
    if(n != MINGrid){
        vec2D bak_f(n+1); for(auto &i:bak_f) i.resize(n+1);
        vec2D f_2(n/2+1); for(auto &i:f_2) i.resize(n/2+1);
        vec2D v_2(n/2+1); for(auto &i:v_2) i.resize(n/2+1);
        // calculate residual
        for(int i=0; i<n; i++){
            bak_f[i][0] = f[i][0] - (1.5*al*n+bl)*v[i][0]+2*al*n*v[i][1]-0.5*al*n*v[i][2];
            bak_f[n][i] = f[n][i] - (1.5*au*n+bu)*v[n][i]+2*au*n*v[n-1][i]-0.5*au*n*v[n-2][i];
        }
        for(int i=1; i<=n; i++){
            bak_f[0][i] = f[0][i] - (1.5*ad*n+bd)*v[0][i]+2*ad*n*v[1][i]-0.5*ad*n*v[2][i];
            bak_f[i][n] = f[i][n] - (1.5*ar*n+br)*v[i][n]+2*ar*n*v[i][n-1]-0.5*ar*n*v[i][n-2];
        }
        for(int i=1; i<n; i++)
            for(int j=1; j<n; j++)
                bak_f[i][j] = f[i][j] -  4*n*n*v[i][j]+n*n*(v[i-1][j]+v[i+1][j]+v[i][j-1]+v[i][j+1]);
        // Restriction
        if(Rest == Operators::fullWeighting){
            for(int i=0; i<=n/2; i++){
                f_2[0][i]   = bak_f[0][2*i];
                f_2[n/2][i] = bak_f[n][2*i];
                f_2[i][0]   = bak_f[2*i][0];
                f_2[i][n/2] = bak_f[2*i][n];
            }
            for(int i=1; i<n/2; i++)
                for(int j=1; j<n/2; j++)
                    f_2[i][j] = 0.0625*(bak_f[2*i-1][2*j-1]+bak_f[2*i+1][2*j+1]+
                                bak_f[2*i-1][2*j+1]+bak_f[2*i+1][2*j-1]) + 
                                0.125*(bak_f[2*i][2*j-1]+bak_f[2*i][2*j+1]+
                                bak_f[2*i-1][2*j]+bak_f[2*i+1][2*j]) + 
                                0.25*bak_f[2*i][2*j];
        }else if(Rest == Operators::injection){
            for(int i=0; i<n/2; i++)
                for(int j=0; j<n/2; j++)
                    f_2[i][j] = bak_f[2*i][2*j];
        }
        // Recursion
        v_2 = V_cycle(w, v_2, f_2, nu1, nu2, Rest, Inter);
        // Prolongation
        if(Inter == Operators::linear){
            for(int i=1; i<=n/2; i++){
                for(int j=1; j<=n/2; j++){
                    v[2*i][2*j]     += v_2[i][j];
                    v[2*i-1][2*j]   += 0.5*v_2[i][j] + 0.5*v_2[i-1][j];
                    v[2*i][2*j-1]   += 0.5*v_2[i][j] + 0.5*v_2[i][j-1];
                    v[2*i-1][2*j-1] += 0.25*v_2[i][j] + 0.25*v_2[i-1][j] +
                                       0.25*v_2[i][j-1] + 0.25*v_2[i-1][j-1];
                }
            }
            for(int i=0; i<n/2; i++){
                v[0][2*i]   += v_2[0][i];
                v[0][2*i+1] += 0.5*v_2[0][i] + 0.5*v_2[0][i+1];
            } v[0][n] += v_2[0][n/2];
            for(int i=1; i<=n/2; i++){
                v[2*i][0]   += v_2[i][0];
                v[2*i-1][0] += 0.5*v_2[i][0] + 0.5*v_2[i-1][0];
            }
        }else if(Inter == Operators::quadratic){
            for(int i=1; i<=n/2; i++){
                for(int j=1; j<=n/2; j++){
                    v[2*i][2*j]     += v_2[i][j];
                    v[2*i-1][2*j]   += 0.5*v_2[i][j] + 0.5*v_2[i-1][j];
                    v[2*i][2*j-1]   += 0.5*v_2[i][j] + 0.5*v_2[i][j-1];
                    v[2*i-1][2*j-1] += 0.25*v_2[i][j] + 0.25*v_2[i-1][j] +
                                       0.25*v_2[i][j-1] + 0.25*v_2[i-1][j-1];
                }
            }
            for(int i=0; i<n/2; i++){
                v[0][2*i]   += v_2[0][i];
                v[0][2*i+1] += 0.5*v_2[0][i] + 0.5*v_2[0][i+1];
            } v[0][n] += v_2[0][n/2];
            for(int i=1; i<=n/2; i++){
                v[2*i][0]   += v_2[i][0];
                v[2*i-1][0] += 0.5*v_2[i][0] + 0.5*v_2[i-1][0];
            }
        }
    }
    // Relax nu2 times
    for(int times=0; times<nu2; times++){
        for(int i=0; i<n; i++){
            temp[i][0] = v[i][0]*(1-w) + v[i][1]*4*al*n*w/(3*al*n+2*bl) - 
                         v[i][2]*al*n*w/(3*al*n+2*bl) + f[i][0]*2*w/(3*al*n+2*bl);
            temp[n][i] = v[n][i]*(1-w) + v[n-1][i]*4*au*n*w/(3*au*n+2*bu) - 
                         v[n-2][i]*au*n*w/(3*au*n+2*bu) + f[n][i]*2*w/(3*au*n+2*bu);
        }
        for(int i=1; i<=n; i++){
            temp[0][i] = v[0][i]*(1-w) + v[1][i]*4*ad*n*w/(3*ad*n+2*bd) -
                         v[2][i]*ad*n*w/(3*ad*n+2*bd) + f[0][i]*2*w/(3*ad*n+2*bd);
            temp[i][n] = v[i][n]*(1-w) + v[i][n-1]*4*ar*n*w/(3*ar*n+2*br) -
                         v[i][n-2]*ar*n*w/(3*ar*n+2*br) + f[i][n]*2*w/(3*ar*n+2*br);
        }
        for(int i=1; i<n; i++)
            for(int j=1; j<n; j++)
                temp[i][j] = v[i][j]*(1-w) + 0.25*w*v[i-1][j] + 0.25*w*v[i+1][j] +
                             0.25*w*v[i][j-1] + 0.25*w*v[i][j+1] + 0.25*h*h*w*f[i][j];
        v = temp;
    }

    return v;
}

vec2D Grid<2>::FMG(real w, vec2D f, int nu1, int nu2, 
            Operators::Restriction Rest, Operators::Interpolation Inter){
    int n = f.size()-1;
    real h = 1.0/n;
    assert((n > 0) && (n&(n-1))==0 && n>=MINGrid);
    // Judgement
    if(n == MINGrid){
        vec2D v(n+1); for(auto &i:v) i.resize(n+1);
        return V_cycle(w,std::move(v),std::move(f),nu1,nu2,Rest,Inter);
    }
    // Restriction and Recursion
    vec2D f_2(n/2+1); for(auto &i:f_2) i.resize(n/+1);
    if(Rest == Operators::fullWeighting){
        for(int i=0; i<=n/2; i++){
            f_2[0][i]   = f[0][2*i];
            f_2[n/2][i] = f[n][2*i];
            f_2[i][0]   = f[2*i][0];
            f_2[i][n/2] = f[2*i][n];
        }
        for(int i=1; i<n/2; i++)
            for(int j=1; j<n/2; j++)
                f_2[i][j] = 0.0625*(f[2*i-1][2*j-1]+f[2*i+1][2*j+1]+
                            f[2*i-1][2*j+1]+f[2*i+1][2*j-1]) + 
                            0.125*(f[2*i][2*j-1]+f[2*i][2*j+1]+
                            f[2*i-1][2*j]+f[2*i+1][2*j]) + 
                            0.25*f[2*i][2*j];
    }else if(Rest == Operators::injection){
        for(int i=0; i<n/2; i++)
            for(int j=0; j<n/2; j++)
                f_2[i][j] = f[2*i][2*j];
    }
    auto v_2 = FMG(w,std::move(f_2),nu1,nu2,Rest,Inter);
    // prolongation
    vec2D v(n+1); for(auto &i:v) i.resize(n+1);
    if(Inter == Operators::linear){
        for(int i=1; i<=n/2; i++){
            for(int j=1; j<=n/2; j++){
                v[2*i][2*j]     = v_2[i][j];
                v[2*i-1][2*j]   = 0.5*v_2[i][j] + 0.5*v_2[i-1][j];
                v[2*i][2*j-1]   = 0.5*v_2[i][j] + 0.5*v_2[i][j-1];
                v[2*i-1][2*j-1] = 0.25*v_2[i][j] + 0.25*v_2[i-1][j] +
                                  0.25*v_2[i][j-1] + 0.25*v_2[i-1][j-1];
            }
        }
        for(int i=0; i<n/2; i++){
            v[0][2*i]   = v_2[0][i];
            v[0][2*i+1] = 0.5*v_2[0][i] + 0.5*v_2[0][i+1];
        } v[0][n] += v_2[0][n/2];
        for(int i=1; i<=n/2; i++){
            v[2*i][0]   = v_2[i][0];
            v[2*i-1][0] = 0.5*v_2[i][0] + 0.5*v_2[i-1][0];
        }
    }else if(Inter == Operators::quadratic){
        for(int i=1; i<=n/2; i++){
            for(int j=1; j<=n/2; j++){
                v[2*i][2*j]     = v_2[i][j];
                v[2*i-1][2*j]   = 0.5*v_2[i][j] + 0.5*v_2[i-1][j];
                v[2*i][2*j-1]   = 0.5*v_2[i][j] + 0.5*v_2[i][j-1];
                v[2*i-1][2*j-1] = 0.25*v_2[i][j] + 0.25*v_2[i-1][j] +
                                  0.25*v_2[i][j-1] + 0.25*v_2[i-1][j-1];
            }
        }
        for(int i=0; i<n/2; i++){
            v[0][2*i]   = v_2[0][i];
            v[0][2*i+1] = 0.5*v_2[0][i] + 0.5*v_2[0][i+1];
        } v[0][n] += v_2[0][n/2];
        for(int i=1; i<=n/2; i++){
            v[2*i][0]   = v_2[i][0];
            v[2*i-1][0] = 0.5*v_2[i][0] + 0.5*v_2[i-1][0];
        }
    }
    // v-cycle
    return V_cycle(w,std::move(v),std::move(f),nu1,nu2,Rest,Inter);
}

void Grid<2>::Poisson_BVP_Multigrid(MultigridMethod method, real w, int nu1, int nu2, 
                           Operators::Restriction Rest, Operators::Interpolation Inter,
                           int max_iter, real epsilon){
#define RESIDUAL                                                                     \
    for(int i=0; i<n; i++){                                                          \
        err[i][0] = abs(RightHandSide_[i][0] -                                       \
                    (1.5*al*n+bl)*v[i][0]+2*al*n*v[i][1]-0.5*al*n*v[i][2]);          \
        err[n][i] = abs(RightHandSide_[n][i] -                                       \
                    (1.5*au*n+bu)*v[n][i]+2*au*n*v[n-1][i]-0.5*au*n*v[n-2][i]);      \
    }                                                                                \
    for(int i=1; i<=n; i++){                                                         \
        err[0][i] = abs(RightHandSide_[0][i] -                                       \
                    (1.5*ad*n+bd)*v[0][i]+2*ad*n*v[1][i]-0.5*ad*n*v[2][i]);          \
        err[i][n] = abs(RightHandSide_[i][n] -                                       \
                    (1.5*ar*n+br)*v[i][n]+2*ar*n*v[i][n-1]-0.5*ar*n*v[i][n-2]);      \
    }                                                                                \
    for(int i=1; i<n; i++)                                                           \
        for(int j=1; j<n; j++)                                                       \
            err[i][j] = abs(RightHandSide_[i][j] -                                   \
                        4*n*n*v[i][j]+n*n*(v[i-1][j]+v[i+1][j]+v[i][j-1]+v[i][j+1]));\
    auto maxlist = std::max_element(err.begin(), err.end(),                          \
                    [](const std::vector<real> &lh, const std::vector<real> &rh){    \
                        return *std::max_element(lh.begin(),lh.end()) <              \
                                *std::max_element(rh.begin(),rh.end());              \
                    });                                                              \
    error = *std::max_element(maxlist->begin(), maxlist->end());

    int n = size_;
    real error;
    auto al = BoundCondi_.leftMixed.first;  auto bl = BoundCondi_.leftMixed.second;
    auto ar = BoundCondi_.rightMixed.first; auto br = BoundCondi_.rightMixed.second;
    auto au = BoundCondi_.upperMixed.first; auto bu = BoundCondi_.upperMixed.second;
    auto ad = BoundCondi_.downMixed.first;  auto bd = BoundCondi_.downMixed.second;
    vec2D v(n+1);   for(auto &i:v)   i.resize(n+1);
    vec2D err(n+1); for(auto &i:err) i.resize(n+1);
    if(method == MultigridMethod::Vcycle){
        for(int i=0; i<max_iter; i++){
            // V-cycle
            v = std::move(V_cycle(w,std::move(v),RightHandSide_,nu1,nu2,Rest,Inter));
            RESIDUAL
            std::cout << "V-cycle执行第" << i+1 << "次，残差为" << error << std::endl;
            if(error < epsilon){break;}
        }

    }else if(method == MultigridMethod::FMG){
        v = std::move(FMG(w,RightHandSide_,nu1,nu2,Rest,Inter));
        RESIDUAL
        std::cout << "FMG执行1次，残差为" << error << std::endl;
    }
    points_ = std::move(v);
    
#undef RESIDUAL
}

void Grid<2>::error_convergece(std::function<real(real,real)> f)
{
    int n = size_;
    real error;
    int max_iter = 100; real epsilon = 1e-8;
    auto al = BoundCondi_.leftMixed.first;  auto bl = BoundCondi_.leftMixed.second;
    auto ar = BoundCondi_.rightMixed.first; auto br = BoundCondi_.rightMixed.second;
    auto au = BoundCondi_.upperMixed.first; auto bu = BoundCondi_.upperMixed.second;
    auto ad = BoundCondi_.downMixed.first;  auto bd = BoundCondi_.downMixed.second;
    vec2D v(n+1);   for(auto &i:v)   i.resize(n+1);
    vec2D err(n+1); for(auto &i:err) i.resize(n+1);
    vec2D solu(n+1); for(auto &i:solu) i.resize(n+1);
    for(int i=0; i<=size_; i++)
        for(int j=0; j<=size_; j++)
            solu[i][j] = f(j*1.0/n, i*1.0/n);
    for(int i=0; i<max_iter; i++){
        // V-cycle
        v = std::move(V_cycle(2.0/3,std::move(v),RightHandSide_,3,3,
                    Operators::fullWeighting,Operators::linear));
        for(int i=0; i<n; i++){                                                          
            err[i][0] = abs(RightHandSide_[i][0] -                                       
                        (1.5*al*n+bl)*v[i][0]+2*al*n*v[i][1]-0.5*al*n*v[i][2]);          
            err[n][i] = abs(RightHandSide_[n][i] -                                       
                        (1.5*au*n+bu)*v[n][i]+2*au*n*v[n-1][i]-0.5*au*n*v[n-2][i]);      
        }                                                                                
        for(int i=1; i<=n; i++){                                                         
            err[0][i] = abs(RightHandSide_[0][i] -                                       
                        (1.5*ad*n+bd)*v[0][i]+2*ad*n*v[1][i]-0.5*ad*n*v[2][i]);          
            err[i][n] = abs(RightHandSide_[i][n] -                                      
                        (1.5*ar*n+br)*v[i][n]+2*ar*n*v[i][n-1]-0.5*ar*n*v[i][n-2]);      
        }                                                                               
        for(int i=1; i<n; i++)                                                           
            for(int j=1; j<n; j++)                                                       
                err[i][j] = abs(RightHandSide_[i][j] -                                   
                            4*n*n*v[i][j]+n*n*(v[i-1][j]+v[i+1][j]+v[i][j-1]+v[i][j+1]));
        auto maxlist = std::max_element(err.begin(), err.end(),                          
                        [](const std::vector<real> &lh, const std::vector<real> &rh){    
                            return *std::max_element(lh.begin(),lh.end()) <              
                                    *std::max_element(rh.begin(),rh.end());              
                        });                                                              
        error = *std::max_element(maxlist->begin(), maxlist->end());
        std::cout << "V-cycle执行第" << i+1 << "次，残差为" << error;

        for(int i=0; i<=n; i++)
            for(int j=0; j<=n; j++)
                err[i][j] = abs(solu[i][j]-v[i][j]);
        auto maxlist2 = std::max_element(err.begin(), err.end(),                          
                        [](const std::vector<real> &lh, const std::vector<real> &rh){    
                            return *std::max_element(lh.begin(),lh.end()) <              
                                    *std::max_element(rh.begin(),rh.end());              
                        });
        std::cout << "; 与解析解的误差为" << *std::max_element(maxlist2->begin(), maxlist2->end())
                  << std::endl;

        if(error < epsilon){break;}
    }


}

#endif