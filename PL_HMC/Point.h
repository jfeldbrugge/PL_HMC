template<size_t dimensions>
class point
{
private:
    std::array<std::complex<double>, dimensions> coords_;
    
public:
    point()
    {
        for(int i = 0; i < dimensions; i++)
        {
            coords_[i] = 0;
        }
    }
    
    point(std::array<std::complex<double>, dimensions> c) : coords_(c)
    {
    }
    
    point(std::initializer_list<std::complex<double>> list)
    {
        size_t n = std::min(dimensions, list.size());
        std::copy_n(list.begin(), n, coords_.begin());
    }
    
    std::complex<double> get(size_t index) const
    {
        return coords_[index];
    }
    
    void assign(std::complex<double> v, size_t index) {
        coords_[index] = v;
    }
    
    point& operator=(const point& a)
    {
        for(size_t i = 0; i < dimensions; ++i)
        {
            coords_[i] = a.get(i);
        }
        return *this;
    }

    point operator+(const point& a) const
    {
        std::array<std::complex<double>, dimensions> tmp;
        for(size_t i = 0; i < dimensions; ++i)
        {
            tmp[i] = coords_[i] + a.get(i);
        }
        return point(tmp);
    }
    
    point operator-(const point& a) const
    {
        std::array<std::complex<double>, dimensions> tmp;
        for(size_t i = 0; i < dimensions; ++i)
        {
            tmp[i] = coords_[i] - a.get(i);
        }
        return point(tmp);
    }

    point operator*(const double &a) const
    {
        std::array<std::complex<double>, dimensions> tmp;
        for(size_t i = 0; i < dimensions; ++i)
        {
            tmp[i] = a * coords_[i];
        }
        return point(tmp);
    }
    
    point operator*(const std::complex<double> &a) const
    {
        std::array<std::complex<double>, dimensions> tmp;
        for(size_t i = 0; i < dimensions; ++i)
        {
            tmp[i] = a * coords_[i];
        }
        return point(tmp);
    }
    
    point operator/(const std::complex<double> &a) const
    {
        std::array<std::complex<double>, dimensions> tmp;
        for(size_t i = 0; i < dimensions; ++i)
        {
            tmp[i] = coords_[i] / a;
        }
        return point(tmp);
    }
    
    bool operator==(const point& a) const
    {
        bool equal = true;
        for (size_t i = 0; i < dimensions; ++i)
        {
            equal = equal && get(i) == a.get(i);
        }
        return equal;
    }
    
    point conjugate() {
        std::array<std::complex<double>, dimensions> tmp;
        for (size_t i = 0; i < dimensions; ++i)
        {
            tmp[i] = std::conj(coords_[i]);
        }
        return point(tmp);
    }
    
    point re() {
        std::array<std::complex<double>, dimensions> tmp;
        for (size_t i = 0; i < dimensions; ++i)
        {
            tmp[i] = std::real(coords_[i]);
        }
        return point(tmp);
    }
    
    point im() {
        std::array<std::complex<double>, dimensions> tmp;
        for (size_t i = 0; i < dimensions; ++i)
        {
            tmp[i] = std::imag(coords_[i]);
        }
        return point(tmp);
    }
};

template<size_t dimensions>
std::ostream& operator<<(std::ostream& out, const point<dimensions>& pt)
{
    out << '(';
    for (size_t i = 0; i < dimensions; ++i)
    {
        if (i > 0)
            out << ",";
        out << pt.get(i);
    }
    out << ')';
    return out;
}

template<size_t dimensions>
double norm(const point<dimensions> &pt)
{
    double norm = 0;
    for(size_t i = 0; i < dimensions; ++i)
    {
        norm = norm + std::real(pt.get(i) * std::conj(pt.get(i)));
    }
    return sqrt(norm);
}

