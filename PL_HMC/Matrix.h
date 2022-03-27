template<size_t dimensions>
class matrix
{
private:
    std::array<std::complex<double>, dimensions * dimensions> elements;
public:
    matrix()
    {
        for(int i = 0; i < dimensions * dimensions; i++)
        {
            elements[i] = 0;
        }
    }
    
    matrix(std::array<std::complex<double>, dimensions * dimensions> c) : elements(c)
    {
    }
    
    matrix(std::initializer_list<std::complex<double>> list)
    {
        size_t n = std::min(dimensions * dimensions, list.size());
        std::copy_n(list.begin(), n, elements.begin());
    }
    
    matrix identity() {
        for(size_t i = 0; i < dimensions; ++i)
        {
            for(size_t j = 0; j < dimensions; ++j)
            {
                if(i == j) {
                    elements[i * dimensions + j] = 1.;
                } else {
                    elements[i * dimensions + j] = 0.;
                }
            }
        }
        return *this;
    }
    
    std::complex<double> get(size_t index1, size_t index2) const
    {
        return elements[index1 * dimensions + index2];
    }
    
    void assign(std::complex<double> v, size_t index1, size_t index2) {
        elements[index1 * dimensions + index2] = v;
    }
    
    matrix& operator=(const matrix& a)
    {
        for(size_t i = 0; i < dimensions; ++i)
        {
            for(size_t j = 0; j < dimensions; ++j)
            {
                elements[i * dimensions + j] = a.get(i, j);
            }
        }
        return *this;
    }
    
    matrix operator+(const matrix& a) const
    {
        std::array<std::complex<double>, dimensions * dimensions> tmp;
        for(size_t i = 0; i < dimensions; ++i)
        {
            for(size_t j = 0; j < dimensions; ++j)
            {
                tmp[i * dimensions + j] = elements[i * dimensions + j] + a.get(i, j);
            }
        }
        return matrix(tmp);
    }
    
    matrix operator-(const matrix& a) const
    {
        std::array<std::complex<double>, dimensions * dimensions> tmp;
        for(size_t i = 0; i < dimensions; ++i)
        {
            for(size_t j = 0; j < dimensions; ++j)
            {
                tmp[i * dimensions + j] = elements[i * dimensions + j] - a.get(i, j);
            }
        }
        return matrix(tmp);
    }
    
    matrix operator*(const double &a) const
    {
        std::array<std::complex<double>, dimensions * dimensions> tmp;
        for(size_t i = 0; i < dimensions * dimensions; ++i)
        {
            tmp[i] = a * elements[i];
        }
        return matrix(tmp);
    }
    
    matrix operator*(const std::complex<double> &a) const
    {
        std::array<std::complex<double>, dimensions * dimensions> tmp;
        for(size_t i = 0; i < dimensions * dimensions; ++i)
        {
            tmp[i] = a * elements[i];
        }
        return matrix(tmp);
    }
    
    matrix operator*(const matrix& M) const
    {
        std::array<std::complex<double>, dimensions * dimensions> tmp;
        for(size_t i = 0; i < dimensions; ++i) {
            for(size_t j = 0; j < dimensions; ++j) {
                std::complex<double> sum = 0;
                for(size_t k = 0; k < dimensions; ++k) {
                    sum = sum + elements[i * dimensions + k] * M.get(k, j);
                }
                tmp[i * dimensions + j] = sum;
            }
        }

        return matrix(tmp);
    }
    
    point<dimensions> operator*(const point<dimensions>& b) const
    {
        std::array<std::complex<double>, dimensions> tmp;
        for(size_t i = 0; i < dimensions; ++i) {
            std::complex<double> sum = 0;
            for(size_t j = 0; j < dimensions; ++j) {
                sum = sum + elements[i * dimensions + j] * b.get(j);
            }
            tmp[i] = sum;
        }
        return point(tmp);
    }
    
    matrix conjugate() {
        std::array<std::complex<double>, dimensions * dimensions> tmp;
        for (int i = 0; i < dimensions * dimensions; i++) {
            tmp[i] = std::conj(elements[i]);
        }
        return matrix(tmp);
    }
    
    std::complex<double> determinant() {
        auto [L, U] = matrix(elements).LU();
        std::complex<double> det = 1.;
        for (int i = 0; i < dimensions; i++) {
            det = det * L.get(i,i) * U.get(i, i);
        }
        return det;
    }
    
    std::tuple<matrix, matrix> LU() const {
        std::array<std::complex<double>, dimensions * dimensions> tmpL;
        std::array<std::complex<double>, dimensions * dimensions> tmpU;
        
        for (int i = 0; i < dimensions; i++) {
            for (int j = 0; j < dimensions; j++) {
                if (j < i)
                    tmpL[j * dimensions + i] = 0;
                else {
                    tmpL[j * dimensions + i] = elements[j * dimensions + i];
                    for (int k = 0; k < i; k++) {
                        tmpL[j * dimensions + i] = tmpL[j * dimensions + i] - tmpL[j * dimensions + k] * tmpU[k * dimensions + i];
                    }
                }
            }
            for (int j = 0; j < dimensions; j++) {
                if (j < i)
                    tmpU[i * dimensions + j] = 0;
                else if (j == i)
                    tmpU[i * dimensions + j] = 1;
                else {
                    tmpU[i * dimensions + j] = elements[i * dimensions + j] / tmpL[i * dimensions + i];
                    for (int k = 0; k < i; k++) {
                        tmpU[i * dimensions + j] = tmpU[i * dimensions + j] - ((tmpL[i * dimensions + k] * tmpU[k * dimensions + j]) / tmpL[i * dimensions  + i]);
                    }
                }
            }
        }
        matrix<dimensions> L(tmpL);
        matrix<dimensions> U(tmpU);
        return {L, U};
    }
    
    matrix transpose() {
        std::array<std::complex<double>, dimensions * dimensions> tmp;
        for (int i = 0; i < dimensions; i++) {
            for (int j = 0; j < dimensions; j++) {
                tmp[i * dimensions + j] = elements[j * dimensions + i];
            }
        }
        return matrix(tmp);
    }
    
    matrix inverse() {
        std::array<std::complex<double>, dimensions * dimensions> tmp;
        matrix A(elements);
        if (A.determinant() != 0.) {
            for (int i = 0; i < dimensions; i++) {
                point<dimensions> b;
                b.assign(1, i);
                point<dimensions> x = solve(A, b);
                for (int j = 0; j < dimensions; j++) {
                    tmp[j * dimensions + i] = x.get(j);
                }
            }
        } else {
            std::cout << "Error: Matrix not invertible!" << std::endl;
        }
        return matrix(tmp);
    }
    
    matrix Cholesky() {
        matrix<dimensions> L, A(elements);
        
        for(int j = 0; j < dimensions; j++) {
            std::complex<double> sum = 0;
            for (int k = 0; k < j; k++) {
                sum += L.get(j, k) * L.get(j,k);
            }
            L.assign(std::sqrt(A.get(j, j) - sum), j, j);

            for (int i = j + 1; i < dimensions; i++) {
                sum = 0;
                for (int k = 0; k < j; k++) {
                    sum += L.get(i, k) * L.get(j, k);
                }
                L.assign((1.0 / L.get(j, j) * (A.get(i, j) - sum)), i, j);
            }
        }
        return L;
    }
};

template<size_t dimensions>
std::ostream& operator<<(std::ostream& out, const matrix<dimensions>& m)
{
    out << '[';
    for (size_t i = 0; i < dimensions; ++i)
    {
        out << '(';
        for (size_t j = 0; j < dimensions; ++j)
        {
            if (j > 0)
                out << ",";
            out << m.get(i, j);
        }
        out << ')';
        if (i != dimensions - 1)
            out << "; ";
    }
    out << ']';
    return out;
}

template<size_t dimensions>
point<dimensions> multiply(const point<dimensions> v, const matrix<dimensions> M)
{
    std::array<std::complex<double>, dimensions> tmp;
    for (int i = 0; i < dimensions; i++) {
        std::complex<double> sum = 0;
        for (int j = 0; j < dimensions; j++) {
            sum = sum + v.get(j) * M.get(j, i);
        }
        tmp[i] = sum;
    }
    return point<dimensions>(tmp);
}

template<size_t dimensions>
matrix<dimensions> identity() {
    std::array<std::complex<double>, dimensions * dimensions> elements;
    for(size_t i = 0; i < dimensions; ++i)
    {
        for(size_t j = 0; j < dimensions; ++j)
        {
            if(i == j) {
                elements[i * dimensions + j] = 1.;
            } else {
                elements[i * dimensions + j] = 0.;
            }
        }
    }
    return matrix<dimensions>(elements);
}

template<size_t dimensions>
point<dimensions> solve(const matrix<dimensions> A, const point<dimensions> b) {
    auto [L, U] = A.LU();
    
    std::array<std::complex<double>, dimensions> tmpY;
    for (int i = 0; i < dimensions; i++) {
        std::complex<double> sum = 0.;
        for (int j = 0; j < dimensions; j++) {
            sum = sum + L.get(i, j) * tmpY[j];
        }
        tmpY[i] = -(sum - b.get(i)) / L.get(i, i);
    }
    point<dimensions> y(tmpY);
    
    std::array<std::complex<double>, dimensions> tmpX;
    for (int i = dimensions - 1; i >= 0; i--) {
        std::complex<double> sum = 0.;
        for (int j = 0; j < dimensions; j++) {
            sum = sum + U.get(i, j) * tmpX[j];
        }
        tmpX[i] = -(sum - y.get(i)) / U.get(i, i);
    }
    point<dimensions> x(tmpX);
    
    return x;
}
