namespace parser {
namespace csv {
void read_iris(std::string fn, std::vector<vec_t> &data, std::vector<int> &label) {
    std::string   s;
    std::ifstream ifs(fn);
    while (std::getline(ifs, s)) {
        vec_t              v;
        std::string        t;
        std::istringstream iss(s);

        // read col1
        std::getline(iss, t, ',');
        v.push_back(std::stof(t));

        // read col2
        std::getline(iss, t, ',');
        v.push_back(std::stof(t));

        // read col3
        std::getline(iss, t, ',');
        v.push_back(std::stof(t));

        // read col4
        std::getline(iss, t, ',');
        v.push_back(std::stof(t));

        data.push_back(v);

        // read label
        std::getline(iss, t, ',');
        label.push_back(std::stoi(t));
    }
}
} // namespace csv
} // namespace parser
