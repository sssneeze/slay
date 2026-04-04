#include "methods/methods_runner.hpp"
#include "methods/symmetric_acceleration_test.hpp"

int main() {
    slau::MethodsRunner::run();
    
    std::cout << "\n\n";
    

    run_symmetric_acceleration_test();
    
    return 0;
}