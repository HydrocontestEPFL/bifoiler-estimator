#include <casadi/casadi.hpp>
#include <iostream>

using namespace casadi;

int main()
{
    SX x = SX::sym("x",1);
    std::cout << "x " << x << std::endl;
}
