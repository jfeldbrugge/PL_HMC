#include "Header.h"

template <size_t dimensions>
std::complex<double> S(const point<dimensions> p) {
    if (dimensions == 1) {
        return (p.get(0) - 1.) * (p.get(0) - 1.) * (p.get(0) + 1.) * (p.get(0) + 1.);
    } else {
        return (p.get(0) - 1.) * (p.get(0) - 1.) * (p.get(0) + 1.) * (p.get(0) + 1.) + p.get(1) * p.get(1);
    }
}

template <size_t dimensions>
std::complex<double> O(const point<dimensions> p) {
    return p.get(0) * p.get(0);
}

template <size_t dimensions>
std::complex<double> O1(const point<dimensions> p) {
    return p.get(0) * p.get(0) * p.get(0) * p.get(0);
}

int main(int argc, const char * argv[]) {
    
    //
    // Parameters
    //
    
    const int dimensions = 1;
    const int N_tau = 10;
    const double tau = 0.1;
    const double Delta_s = 0.02;
    const int N_s = 200;
    const int N_samples = 1000000;
    const double m = 2.;

    matrix<dimensions> M = identity<dimensions>() * m;
    matrix<dimensions> M_inv = M.inverse(), M_Cholesky = M.Cholesky();
    
    std::srand(0);
    std::default_random_engine gen;
    std::normal_distribution<double> normal(0., 1);
    
    std::cout << "Picard-Lefschetz Hamiltonian Monte Carlo with " << N_samples << " samples" << std::endl;
    
    //
    // Sample the points
    //
    
    std::vector<point<dimensions>> xi(N_samples);
    std::cout << "Progress:"; std::cout.flush();
    point<dimensions> x;
    int count = 0;
    int reject = 0;
    while (count < N_samples) {
        point<dimensions> p;
        for (int i = 0; i < dimensions; i++) p.assign(normal(gen), i);
        p = M_Cholesky * p;

        double H0 = hamiltonian(x, p, M_inv, tau, N_tau);
        auto [X, H1] = leapfrog(x, p, Delta_s, N_s, tau, N_tau, M_inv);
        double delta_H = H1 - H0;

        if(!std::isnan(std::real(X.get(0))) && (std::real(X * X.conjugate()) < 1000) && (uniform() < std::min(1., std::exp(- delta_H)))) {
            xi[count] = X;
            x = X;
            count++;
            if ((count % (N_samples / 10) == 0)) {
                std::cout << "."; std::cout.flush();
            }
        } else {
            reject++;
        }
    }
    std::cout << std::endl;
    std::cout << "Acceptance/Rejection rate: " << double(N_samples) / double(N_samples + reject) << "/" <<
                                                  double(reject) / double(N_samples + reject)    << std::endl;
    writeB(xi, "xi_" + std::to_string(dimensions) + ".bin");
    
    //
    // Evaluate the expectation value
    //
    
    std::complex<double> mean = expectation(&O, xi, tau, N_tau);
    std::complex<double> var = variance(&O, mean, xi, tau, N_tau);

    std::cout << "E[O(x)] = " << mean << " +/- " << std::real(std::sqrt(var)) / std::sqrt(double(N_samples)) << std::endl;
    
    std::cout << "Done!" << std::endl;
    return 0;
}
