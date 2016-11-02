#include "./util.hpp"

using namespace cc;

#include "./parser.hpp"
#include "./pca.hpp"

int main(int argc, char *argv[]) {
    try {
        // アイリスデータを読み込む
        std::vector<vec_t> data;
        std::vector<int>   label;
        parser::csv::read_iris("./iris.csv", data, label);

        // 主成分分析する
        std::vector<vec_t> flatten;
        matrix::pca(data, flatten);

        // プロット用に出力する
        std::ofstream ofs("output");
        for (size_t i = 0; i < flatten.size(); i++) {
            ofs << format_str("%f %f %d\n", flatten[i][0], flatten[i][1], label[i]);
        }

    } catch (const std::exception &e) {
        std::cerr << colorant('y', format_str("error: %s", e.what())) << std::endl;
    }
}


