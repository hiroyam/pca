namespace matrix {
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

void pca(std::vector<vec_t> &data, std::vector<vec_t> &flatten) {
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
	matrix::jacobi(covar, e, ev);

	// 1番目に大きい固有値のインデクスを調べる
	int i1 = max_index(e);

	// 再び選ばれないように0で埋めておく
	e[i1] = 0.0f;

	// 2番めに大きい固有値のインデクスを調べる
	int i2 = max_index(e);

	// 選んだ固有値の固有ベクトルで潰す
	flatten.clear();
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
}

} // matrix
