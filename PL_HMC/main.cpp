#include "Header.h"

template <size_t dimensions>
std::complex<double> S(const point<dimensions> p) {
    return (p.get(0) - 1.) * (p.get(0) - 1.) * (p.get(0) + 1.) * (p.get(0) + 1.);
}

template <size_t dimensions>
std::complex<double> O(const point<dimensions> p) {
    return p.get(0) * p.get(0);
}
int main(int argc, const char * argv[]) {
    // Parameters
    const int dimensions = 1;
    const int N_tau = 10;
    const double tau = 0.1;
    const double Delta_s = 0.02;
    const int N_s = 200;
    const int N_samples = 1000000;
    const double m = 1.;
    
    matrix<3> M({10., 1., 0., 2., 1., 3., 3., 1., 1.});
    
    std::srand(0);
    std::default_random_engine gen;
    std::normal_distribution<double> normal(0., m);
    
    //
    // Sample the points
    //
    
    std::cout << "Picard-Lefschetz Hamiltonian Monte Carlo with " << N_samples << " samples" << std::endl;
    std::vector<point<dimensions>> xi(N_samples);
    std::cout << "Progress:";
    point<dimensions> x({0.1});
    int count = 0;
    int reject = 0;
    while (count < N_samples) {
        point<dimensions> p({normal(gen)});
        
        double H0 = Hamiltonian(x, p, m, tau, N_tau);
        auto [X, H1] = leapfrog(x, p, Delta_s, N_s, tau, N_tau, m);
        double delta_H = H1 - H0;
        if(uniform() < std::min(1., std::exp(- delta_H))) {
            xi[count] = X;
            x = X;
            count++;
            if ((count % (N_samples / 100) == 0)) {
                std::cout << "."; std::cout.flush();
            }
        } else {
            reject++;
        }
    }
    std::cout << std::endl;
    std::cout << "Rejection rate: " << double(reject) / double(N_samples) * 100. << "%" << std::endl;
    writeB(xi, "xi.bin");
    
    //
    // Evaluate the expectation value
    //
    
    std::complex<double> numerator = 0., denominator = 0.;
    for (int i = 0; i < xi.size(); i++) {
        point<dimensions> z, f;
        matrix<dimensions> J;
        flow(xi[i], z, J, f, tau, N_tau);
        numerator = numerator + O(z) * std::exp(I * H(z)) * determinant(J);
        denominator = denominator + std::exp(I * H(z)) * determinant(J);
    }
    std::cout << "E[O(x)] = " << numerator / denominator << std::endl;
    
    std::cout << "Done!" << std::endl;
    return 0;
}
