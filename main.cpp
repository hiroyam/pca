#include "./util.hpp"

using namespace cc;

#include "./matrix.hpp"
#include "./parser.hpp"

int main(int argc, char *argv[]) {
    try {
        std::vector<vec_t> data;
        std::vector<int>   label;

        // アイリスデータを読み込む
        parser::csv::read_iris("./iris.csv", data, label);

        int N   = data.size();
        int dim = data[0].size();

        // 平均を計算する
        vec_t mean(dim);
        for (auto v : data) {
            for (int i = 0; i < dim; i++) {
                mean[i] += v[i] / N;
            }
        }

        // データの平均を0にする
        for (auto &v : data) {
            for (int i = 0; i < dim; i++) {
                v[i] -= mean[i];
            }
        }

        // 分散共分散行列を計算する
        vec_t covar(dim * dim);
        for (auto v : data) {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    covar[i + j * dim] += v[i] * v[j] / N;
                }
            }
        }

        // 固有値を計算する
        vec_t e;
        vec_t ev;
        matrix::eigen::jacobi(covar, e, ev);

        // 1番目に大きい固有値のインデクスを調べる
        int i1 = max_index(e);

        // 再び選ばれないように0で埋めておく
        e[i1] = 0.0f;

        // 2番めに大きい固有値のインデクスを調べる
        int i2 = max_index(e);

        // 選んだ固有値の固有ベクトルで潰す
        std::vector<vec_t> flatten;
        for (auto v : data) {
            float x = 0;
            float y = 0;
            for (int i = 0; i < dim; i++) {
                x += v[i] * ev[i1 + i * dim];
                y += v[i] * ev[i2 + i * dim];
            }
            vec_t f{x, y};
            flatten.push_back(f);
        }

        // プロット用に出力する
        std::ofstream ofs("output");
        for (int i = 0; i < N; i++) {
            ofs << format_str("%f %f %d\n", flatten[i][0], flatten[i][1], label[i]);
        }

    } catch (const std::exception &e) {
        std::cerr << colorant('y', format_str("error: %s", e.what())) << std::endl;
    }
}

