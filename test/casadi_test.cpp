#include <casadi/casadi.hpp>
#include <math.h>
#include <iostream>

using namespace casadi;

int main()
{
    SX x = SX::sym("x",1);
    std::cout << "x " << x << std::endl;
    auto f = SX::vertcat({sin(x), pow(x, 2)});
    std::cout << "f " << f << std::endl;
    std::cout << "Jacobian(f) " << SX::jacobian(f,x) << std::endl;
}
