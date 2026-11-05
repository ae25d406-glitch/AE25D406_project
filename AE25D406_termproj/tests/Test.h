#ifndef TEST_H
#define TEST_H

#include <cmath>
#include <iostream>
#include <string>

class Test
{
private:

    int passed = 0;
    int failed = 0;

public:

    void expectTrue
    (
        bool condition,
        const std::string& message
    )
    {
        if(condition)
        {
            std::cout << "[PASS] "
                      << message
                      << std::endl;

            passed++;
        }
        else
        {
            std::cout << "[FAIL] "
                      << message
                      << std::endl;

            failed++;
        }
    }

    void expectEqual
    (
        double actual,
        double expected,
        double tolerance,
        const std::string& message
    )
    {
        if(std::fabs(actual-expected)<tolerance)
        {
            std::cout << "[PASS] "
                      << message
                      << std::endl;

            passed++;
        }
        else
        {
            std::cout << "[FAIL] "
                      << message
                      << std::endl;

            std::cout
                << " expected : "
                << expected
                << std::endl;

            std::cout
                << " actual   : "
                << actual
                << std::endl;

            failed++;
        }
    }

    void summary()
    {
        std::cout
            << "\n==========================\n";

        std::cout
            << "Passed : "
            << passed
            << std::endl;

        std::cout
            << "Failed : "
            << failed
            << std::endl;

        std::cout
            << "==========================\n";

        if(failed==0)
        {
            std::cout
                << "ALL TESTS PASSED\n";
        }
    }

    bool success() const
    {
        return failed==0;
    }

};

#endif