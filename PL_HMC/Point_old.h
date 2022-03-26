template<typename coordinate_type, size_t dimensions>
class point
{
private:
    std::array<coordinate_type, dimensions> coords_;
    
public:
    point()
    {
        for(int i = 0; i < dimensions; i++)
        {
            coords_[i] = 0;
        }
    }
    
    point(std::array<coordinate_type, dimensions> c) : coords_(c)
    {
    }
    
    point(std::initializer_list<coordinate_type> list)
    {
        size_t n = std::min(dimensions, list.size());
        std::copy_n(list.begin(), n, coords_.begin());
    }
    
    coordinate_type get(size_t index) const
    {
        return coords_[index];
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
        std::array<coordinate_type, dimensions> tmp;
        for(size_t i = 0; i < dimensions; ++i)
        {
            tmp[i] = coords_[i] + a.get(i);
        }
        return point(tmp);
    }
    
    point operator-(const point& a) const
    {
        std::array<coordinate_type, dimensions> tmp;
        for(size_t i = 0; i < dimensions; ++i)
        {
            tmp[i] = coords_[i] - a.get(i);
        }
        return point(tmp);
    }

    point operator*(const coordinate_type &a) const
    {
        std::array<coordinate_type, dimensions> tmp;
        for(size_t i = 0; i < dimensions; ++i)
        {
            tmp[i] = a * coords_[i];
        }
        return point(tmp);
    }
    
    point operator/(const coordinate_type &a) const
    {
        std::array<coordinate_type, dimensions> tmp;
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
};

template<typename coordinate_type, size_t dimensions>
std::ostream& operator<<(std::ostream& out, const point<coordinate_type, dimensions>& pt)
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

template<typename coordinate_type, size_t dimensions>
double norm(const point<coordinate_type, dimensions> &pt)
{
    double norm = 0;
    for(size_t i = 0; i < dimensions; ++i)
    {
        norm = norm + std::real(pt.get(i) * std::conj(pt.get(i)));
    }
    return sqrt(norm);
}

