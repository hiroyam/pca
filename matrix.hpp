namespace matrix {
namespace eigen {
/**
 * @param a  実数対称行列
 * @param e  固有値
 * @param ev 各列が固有ベクトル
 */
void jacobi(const vec_t &a, vec_t &e, vec_t &ev) {
    int N   = a.size();
    int dim = std::sqrt(N);

    if (N != dim * dim) throw std::runtime_error("failed to jacobi: invalid vector");

    // aのサイズに合わせる
    e  = a;
    ev = a;

    // evを単位行列にする
    for (int j = 0; j < dim; j++) {
        for (int i = 0; i < dim; i++) {
            ev[i + j * dim] = (i == j) ? 1.0f : 0.0f;
        }
    }

    int loop_count = 0;
    while (true) {
        // 非対角成分のうち絶対値が最大の要素を探す
        int   p = 0;
        int   q = 0;
        float m = 0;
        for (int r = 0; r < dim; r++) {
            for (int c = 0; c < dim; c++) {
                if (r >= c) continue;

                float val = std::abs(e[c + r * dim]);
                if (val > m) {
                    p = r;
                    q = c;
                    m = val;
                }
            }
        }
        // 絶対値が最大の要素が許容誤差以内になったら終了する
        const float eps = 0.001f;
        if (m <= eps) break;

        // 無限ループに陥ったら終了する
        if (++loop_count > 1000) {
            throw std::runtime_error("failed to calculate eigen");
        }

        float app = e[p + p * dim];
        float aqq = e[q + q * dim];
        float apq = e[q + p * dim];

        float alpha = (app - aqq) / 2.0f;
        float beta  = -apq;
        float gamma = std::abs(alpha) / std::sqrt(alpha * alpha + beta * beta);

        float ct = std::sqrt((1.0f + gamma) / 2.0f);
        float st = std::sqrt((1.0f - gamma) / 2.0f);
        if (alpha * beta < 0) st = -st;

        // p列とq列の更新
        for (int r = 0; r < dim; r++) {
            float tmp_p = ct * e[p + r * dim] - st * e[q + r * dim];
            float tmp_q = st * e[p + r * dim] + ct * e[q + r * dim];
            e[p + r * dim] = tmp_p;
            e[q + r * dim] = tmp_q;
        }

        // 対称にする
        for (int r = 0; r < dim; r++) {
            e[r + p * dim] = e[p + r * dim];
            e[r + q * dim] = e[q + r * dim];
        }

        // pp, qq, pq, qpの更新
        e[p + p * dim] = ct * ct * app + st * st * aqq - 2 * st * ct * apq;
        e[q + q * dim] = st * st * app + ct * ct * aqq + 2 * st * ct * apq;
        e[p + q * dim] = ct * st * (app - aqq) + (ct * ct - st * st) * apq;
        e[q + p * dim] = ct * st * (app - aqq) + (ct * ct - st * st) * apq;

        // 対角化行列の更新
        for (int r = 0; r < dim; r++) {
            float tmp_p = ct * ev[p + r * dim] - st * ev[q + r * dim];
            float tmp_q = st * ev[p + r * dim] + ct * ev[q + r * dim];
            ev[p + r * dim] = tmp_p;
            ev[q + r * dim] = tmp_q;
        }
    }

    // eの対角成分が固有値なのでそこだけ切り出す
    for (int i = 0; i < dim; i++) {
        e[i] = e[i + i * dim];
    }
    e.resize(dim);
}
} // namespace eigen
} // namespace matrix
