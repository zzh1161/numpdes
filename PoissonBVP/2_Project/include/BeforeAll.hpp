#ifndef __PRO2_INFO_HPP__
#define __PRO2_INFO_HPP__

using real = double;

enum class Domin{
    regular = 0, 
    irregular
};

enum class BoundCondiType{
    Dirichlet = 0,
    Neumann,
    Mixed
};

enum class MultigridMethod{
    Vcycle = 0,
    FMG
};

namespace Operators
{
    enum Restriction{
        injection = 0,
        fullWeighting
    };
    enum Interpolation{
        linear = 0,
        quadratic
    };
};

#endif