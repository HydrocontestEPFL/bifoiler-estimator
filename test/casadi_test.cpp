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

    SX q = SX::vertcat({0, 1, 2, 3});
    std::cout << "q(2:4) " << q(Slice(2,4)) << std::endl;
    std::cout << "q(1) " << q(1) << std::endl;

    Function dyn = Function("dynamics", {x}, {f});
    std::cout << dyn << "\n";
    dyn.generate();

    SX A = SX::sym("A", 3, 3);
    SX A_sliced = A(Slice(0,2),Slice(0,2));
    std::cout << "A " << A << "\n";
    std::cout << "A(0:2,0:2) " << A_sliced << "\n";
}
