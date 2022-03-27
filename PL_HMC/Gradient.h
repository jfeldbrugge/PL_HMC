template <size_t dimensions>
point<dimensions> gradS(const point<dimensions> p)
{
    const double eps = 0.000001;
    std::array<std::complex<double>, dimensions> grad;
    
    for (int i = 0; i < dimensions; i++) {
        std::array<std::complex<double>, dimensions> temp;
        temp[i] = eps;
        point<dimensions> Peps(temp);
        grad[i] = (S(p + Peps) - S(p - Peps)) / (2. * eps);
    }
    
    return point<dimensions>(grad);
}

template <size_t dimensions>
matrix<dimensions> hessS(const point<dimensions> p)
{
    const double eps = 0.0001;
    std::array<std::complex<double>, dimensions * dimensions> hessian;

    for (int i = 0; i < dimensions; i++) {
        for (int j = 0; j < dimensions; j++) {
            
            if (i == j) {
                std::array<std::complex<double>, dimensions> temp;
                temp[i] = eps;
                point<dimensions> Peps(temp);
                hessian[i * dimensions + j] = (S(p + Peps) -2. * S(p) + S(p - Peps)) / (eps * eps);
                
            } else {
                std::array<std::complex<double>, dimensions> temp1, temp2;
                temp1[i] = eps; temp2[j] = eps;
                point<dimensions> Peps1(temp1), Peps2(temp2);
                hessian[i * dimensions + j] = (S(p + Peps1 + Peps2) - S(p + Peps1 - Peps2) - S(p - Peps1 + Peps2) + S(p - Peps1 - Peps2)) / (4. * eps * eps);
            }
            
//            std::cout << hessian[i * dimensions + j] << std::endl;
        }
    }

    return matrix<dimensions>(hessian);
}
