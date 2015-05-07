#pragma once

#include <complex>
#include <vector>

class Kuramoto
{
    public:
        explicit Kuramoto(size_t network_size
                         ,double stepsize
                         ,std::vector<double> couplings
                         ,std::vector<double> phases
                         ,std::vector<double> nat_frequencies);

        // Delegated constructor
        explicit Kuramoto(size_t network_size
                         ,double stepsize
                         ,double coupling_strength
                         ,std::vector<double> phases
                         ,std::vector<double> nat_frequencies);

        ~Kuramoto(){}

        // Function describing the state and evolution of the system
        std::complex<double> calculate_order_param();
        std::tuple<double, double> polar(const std::complex<double> &z);
        void make_step();
        std::vector<double> run(int steps);


    private:
        size_t m_size;
        double m_stepsize;
        int m_stepcount {0};
        std::vector<double> m_couplings;
        std::vector<double> m_phases;
        std::vector<double> m_nat_frequencies;
};
