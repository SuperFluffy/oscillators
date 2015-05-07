#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <random>
#include <stdexcept>
#include <vector>

#include "kuramoto.hh"

Kuramoto::Kuramoto(size_t network_size
         ,double stepsize
         ,std::vector<double> couplings
         ,std::vector<double> phases
         ,std::vector<double> nat_frequencies
         )
:   m_size{network_size}, m_stepsize{stepsize}, m_couplings(couplings),
    m_phases(phases), m_nat_frequencies(nat_frequencies)
{
    if (m_couplings.size() != network_size) {
        throw std::length_error("Size of coupling vector not equal to sytem size.");
    }
    if (m_phases.size() != network_size) {
        throw std::length_error("Size of phase vector not equal to sytem size.");
    }
    if (m_nat_frequencies.size() != network_size) {
        throw std::length_error("Size of frequency vector not equal to sytem size.");
    }
}

// Delegated constructor
Kuramoto::Kuramoto(size_t network_size
         ,double stepsize
         ,double coupling_strength
         ,std::vector<double> phases
         ,std::vector<double> nat_frequencies
         )
: Kuramoto{network_size,stepsize
          ,std::vector<double>(network_size,coupling_strength)
          ,phases
          ,nat_frequencies}
{
//
}

void Kuramoto::make_step()
{
    std::vector<double> k1(m_size, 0);

    auto k = k1.begin();
    auto kend = k1.end();
    auto phase = m_phases.begin();

    for (; k != kend; ++k, ++phase)
    {
        auto p = m_phases.begin();
        auto pend = m_phases.end();
        auto c = m_couplings.begin();
        for (; p != pend ; ++p, ++c)
        {
            *k += *c * sin(*p - *phase);
        }
    }

    auto p = m_phases.begin();
    auto pend = m_phases.end();
    auto f = m_nat_frequencies.begin();
    k = k1.begin();
    for (; p != pend ; ++p, ++f, ++k)
    {
        *p += m_stepsize * (*f + *k);
    }

    ++m_stepcount;

    // Normalize the phases to be in the range [0,Pi)
    using namespace std::placeholders;

    //auto mod_fn = std::bind((double(*)(double,double))std::fmod, _1, M_PI);
    //auto mod_fn = std::bind(static_cast<double(*)(double,double)>(std::fmod), _1, M_PI);
    auto mod_fn = [](double d) { return std::fmod(d, M_PI); };
    //auto mod_fn = std::bind<double(double,double)>(&std::fmod, _1, M_PI);
    std::transform(m_phases.begin(), m_phases.end(), m_phases.begin(), mod_fn);
}

std::complex<double> Kuramoto::calculate_order_param()
{
    // ∑ᵢ exp(ⅈ*φᵢ), with φᵢ the phase of the i-th oscillator
    std::complex<double> order = std::accumulate(m_phases.begin(), m_phases.end(), 0.0,
                                         [](double &sum, const double &p)
                                         { return sum + std::exp(p); }
                                         );
    order /= m_size;
    return order;
}

std::tuple<double,double> Kuramoto::polar(const std::complex<double> &z)
{
    return std::tuple<double,double>(std::abs(z), std::arg(z));
}
