template <size_t dimensions>
double hamiltonian(point<dimensions> z, point<dimensions> p, const matrix<dimensions> M_inv) {
    return std::real(M_inv * p * p / 2.) - h(z);
}

template <size_t dimensions>
double hamiltonian(point<dimensions> x, point<dimensions> p, const matrix<dimensions> M_inv, const double tau, const int N_tau) {
    point<dimensions> z, f;
    matrix<dimensions> J;
    flow(x, z, J, f, tau, N_tau);
    return hamiltonian(z, p, M_inv);
}
    
template <size_t dimensions>
std::tuple<point<dimensions>, double> leapfrog(const point<dimensions> x, const point<dimensions> p,
              const double Delta_s, const int N_s, const double tau, const int N_tau, const matrix<dimensions> Minv) {
    point<dimensions> X = x, P = p;
    point<dimensions> F = force(X, tau, N_tau);
    for (int i = 1; i < N_s + 1; i++) {
        // Kick
        P = P + F * Delta_s / 2.;
        // Drift
        X = X + Minv * P * Delta_s;
        //Kick
        F = force(X, tau, N_tau);
        P = P + F * Delta_s / 2.;
    }
    return {X, hamiltonian(X, P, Minv, tau, N_tau)};
}

template <size_t dimensions>
std::complex<double> expectation(std::complex<double> (*OO)(point<dimensions>), const std::vector<point<dimensions>> &xi, const double tau, const int N_tau) {
    std::complex<double> numerator = 0., denominator = 0.;
    for (int i = 0; i < xi.size(); i++) {
        point<dimensions> z, f;
        matrix<dimensions> J;
        flow(xi[i], z, J, f, tau, N_tau);
        numerator   = numerator   + OO(z) * std::exp(I * H(z)) * J.determinant();
        denominator = denominator +         std::exp(I * H(z)) * J.determinant();
    }
    return numerator / denominator;
}
