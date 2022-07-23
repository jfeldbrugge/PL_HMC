template<size_t dimensions>
void writeB(std::vector<point<dimensions>> xi, std::string fileName) {
    std::ofstream file; file.open(fileName, std::ios::binary);
    if(file.is_open()) {
        for(int index = 0; index < xi.size(); index++) {
            for (int i = 0; i < dimensions; i++) {
                double g = std::real(xi[index].get(i));
                file.write((char*) &g, sizeof(double));
            }
        }
    } else {
        std::cout << "Could not open " << fileName << std::endl;
    }
    file.close();
}
