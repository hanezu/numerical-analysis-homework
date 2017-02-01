// u = x^2 + y^2 + 4*t
// f = 0
// u_D = u
// u^0(x,y) = x^2 + y^2
// Region is circle
// by: Huachun Zhu


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

//#define GAMMA (1.0)
#define GAMMA (1/2.0)

// 矩形領域 [X0, X0 + L] \times [Y0, Y0 + L] に格子点を作る
// 格子の基準点
#define X0 (-1.0)
#define Y0 (-1.0)
// 格子の縦および横の長さ
#define L (2.0)
//time length
#define T (2.0)
// L あたりの分割数
int N = 50;
// T あたりの分割数
int r = 5;

// 分割幅
#define h (L / N)
#define tau (T / r)



// ヤコビ法の反復終了の閾値
#define TOLERANCE (1.0e-13)
// ヤコビ法の最大反復回数
#define MAXITER (20000)
// 出力時に間引くデータ間隔
#define CULL (1)

FILE *fd;
FILE *fr;

// change indexes to coordinates
double i2x(int i) {
    return X0 + i * h;
}

double j2y(int j) {
    return Y0 + j * h;
}

double k2t(int k) {
    return k * tau;
}

// 内部節点に連番をふるための配列
int **indices;

// 方程式の右辺
double f(double x, double y, double t) {
    return 0.0;
}

// 境界条件
double Dirichlet_data(double x, double y, double t) {
    return x * x + y * y + 4 * t;
}

double Start_data(double x, double y) {
    return x * x + y * y;
}


// 境界のパラメータ表示
// 凸領域でない場合，上下左右それぞれの表示が一意に定まるとは限らないことに注意
double y_top(double x) {
    return sqrt(1.0 - x * x);
}

double y_bottom(double x) {
    return -sqrt(1.0 - x * x);
}

double x_right(double y) {
    return sqrt(1.0 - y * y);
}

double x_left(double y) {
    return -sqrt(1.0 - y * y);
}

// 内部節点かどうか判定する．
// 境界は false を返す．
bool is_internal(int i, int j) {
    double x = X0 + i * h;
    double y = Y0 + j * h;
    if (y >= y_top(x)) return false;
    if (y <= y_bottom(x)) return false;
    if (x >= x_right(y)) return false;
    if (x <= x_left(y)) return false;
    return true;
}

// 内部節点 P_ij における lambda_X (X=B,C,D,E) の値を返す．
// 凸領域で上下左右の媒介変数表示が与えられているため以下のように
// 簡単に求まるが，一般にこの方法で求められるとは限らない．
// （つまり，より丁寧な場合分けが必要である．）
// 右側
double lambda_B(int i, int j) {
    if (is_internal(i + 1, j)) return 1.0;
    double lambda = x_right(Y0 + j * h) - (X0 + i * h);//this is [0,h) !!! why not divide it by h
    return lambda / h;

}

// 左側
double lambda_C(int i, int j) {
    if (is_internal(i - 1, j)) return 1.0;
    double lambda = (X0 + i * h) - x_left(Y0 + j * h);
    return lambda / h;
}

// 上側
double lambda_D(int i, int j) {
    if (is_internal(i, j + 1)) return 1.0;
    double lambda = y_top(X0 + i * h) - (Y0 + j * h);
    return lambda / h;
}

// 下側
double lambda_E(int i, int j) {
    if (is_internal(i, j - 1)) return 1.0;
    double lambda = (Y0 + j * h) - y_bottom(X0 + i * h);
    return lambda / h;
}

// average of the lambdas is useful
double avg_lambda(int i, int j) {
    if (is_internal(i, j)) {
        return (lambda_B(i, j) + lambda_C(i, j)) * (lambda_D(i, j) + lambda_E(i, j)) / 4.0;
    } else {
        fprintf(stderr, "don't use sum_lambda with external points!");
        exit(EXIT_FAILURE);
    }
}

// この関数は下で定義する．
double *allocate_real_vector(size_t);


void test(void) {

}

void run_main(int my_N, int my_r);

int main(int argc, char *argv[]) {
    bool is_test = false;
    if (is_test) {
        test();
    } else {
//        int N_array[] = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190};
//        int N_r_array[1][2] = {{100, 100}};
//        int len_N_r_array = sizeof(N_r_array) / sizeof(N_r_array[0]);
        fr = fopen("error_fit.csv", "w+");
        fprintf(fr, "N r h tau error\n");

        for (int i = 10; i <= 100; i += 10) {
            for (int j = 10; j <= 100; j += 10) {
                run_main(i, j);
            }
        }
        fclose(fr);
    }

    return 0;
}

void run_main(int my_N, int my_r) {
    N = my_N;
    r = my_r;

    indices = malloc((N + 1) * sizeof(*indices)); //init indices row
    fd = fopen("debug.txt", "w");
    // 内部節点の数 n を数える（未知数の数）
    int n = 0;
    for (int i = 0; i <= N; ++i) {
        indices[i] = malloc((N + 1) * sizeof(*indices[i]));
        for (int j = 0; j <= N; ++j) {
            if (is_internal(i, j)) {
                // 節点に順番に番号を振って格納しておく
                indices[i][j] = n;
                ++n;
            } else {
                indices[i][j] = -1;
            }
        }
    }

    // // // // // // // // // * キ ケ ン * // // // // // // // // // //
    // 事前に未知数の数 n が分からないので配列は動的確保する．
    // 動的確保した配列は free() を使って解放するのを忘れないこと．
    double **uh = malloc((r + 1) * sizeof(*uh));
    double *uh_old = allocate_real_vector(n); // notice that it initialize anyway whenever k move on
//    double *uh_new = allocate_real_vector(n);
    double *uh_new;
    double *Fh = allocate_real_vector(n);

    // init u[t=0]
    uh[0] = allocate_real_vector(n);
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            if (is_internal(i, j)) {
                uh[0][indices[i][j]] = Start_data(i2x(i), j2y(j));
            }
        }
    }

    // we will solve N * N equations for every time k
    // and after solving, we move the k forward
    for (int k = 0; k <= r - 1; k++) {
        uh_new = allocate_real_vector(n);

        // 右辺ベクトルを初期化する
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                if (is_internal(i, j)) {
                    int m = indices[i][j];
                    double x = X0 + i * h;
                    double y = Y0 + j * h;
                    double t = k2t(k);
                    // 正則な節点の場合，C = 1 となることに注意
                    double L_B = 1.0 / lambda_B(i, j);
                    double L_C = 1.0 / lambda_C(i, j);
                    double L_D = 1.0 / lambda_D(i, j);
                    double L_E = 1.0 / lambda_E(i, j);
                    double C = (lambda_B(i, j) + lambda_C(i, j)) * (lambda_D(i, j) + lambda_E(i, j)) / 4.0;
                    Fh[m] = C * h * h * (f(x, y, t + tau) * GAMMA + f(x, y, t) * (1 - GAMMA)) // part of f
                            + (C * h * h / tau - (1 - GAMMA) * (L_B + L_C + L_D + L_E)) * uh[k][m];
                    // 境界条件から既知の数を右辺ベクトルに足し上げる．
                    if (!is_internal(i + 1, j)) { // 右側
                        Fh[m] += ((1 - GAMMA) * Dirichlet_data(x_right(y), y, t) +
                                  GAMMA * Dirichlet_data(x_right(y), y, t + tau)) / lambda_B(i, j);
                    } else { //use u[k-1] data
                        Fh[m] += (1 - GAMMA) * uh[k][indices[i + 1][j]] / lambda_B(i, j);
                    }
                    if (!is_internal(i - 1, j)) { // 左側
                        Fh[m] += ((1 - GAMMA) * Dirichlet_data(x_left(y), y, t) +
                                  GAMMA * Dirichlet_data(x_left(y), y, t + tau)) / lambda_C(i, j);
                    } else { //use u[k-1] data
                        Fh[m] += (1 - GAMMA) * uh[k][indices[i - 1][j]] / lambda_C(i, j);
                    }
                    if (!is_internal(i, j + 1)) { // 上側
                        Fh[m] += ((1 - GAMMA) * Dirichlet_data(x, y_top(x), t) +
                                  GAMMA * Dirichlet_data(x, y_top(x), t + tau)) / lambda_D(i, j);
                    } else { //use u[k-1] data
                        Fh[m] += (1 - GAMMA) * uh[k][indices[i][j + 1]] / lambda_D(i, j);
                    }
                    if (!is_internal(i, j - 1)) { // 下側
                        Fh[m] += ((1 - GAMMA) * Dirichlet_data(x, y_bottom(x), t) +
                                  GAMMA * Dirichlet_data(x, y_bottom(x), t + tau)) / lambda_E(i, j);
                    } else { //use u[k-1] data
                        Fh[m] += (1 - GAMMA) * uh[k][indices[i][j - 1]] / lambda_E(i, j);
                    }

                }
            }
        }

        // ヤコビ法による反復計算で方程式 Ah * uh = Fh を解く
        // 反復の初期値は D^{-1} * Fh に近いものでよい
        for (int m = 0; m < n; ++m) {
            uh_old[m] = Fh[m] * h * h / 4.0;
        }

        for (int l = 0; l < MAXITER; ++l) {
            // 疎行列 H = -D^{-1} * R は非零要素を非対角に高々 4 つもつ
            for (int i = 0; i <= N; ++i) {
                for (int j = 0; j <= N; ++j) {
                    if (is_internal(i, j)) {
                        int m = indices[i][j];
                        double L_B = 1.0 / lambda_B(i, j);
                        double L_C = 1.0 / lambda_C(i, j);
                        double L_D = 1.0 / lambda_D(i, j);
                        double L_E = 1.0 / lambda_E(i, j);
                        double sum_L = L_B + L_C + L_D + L_E;
                        double C = (lambda_B(i, j) + lambda_C(i, j)) * (lambda_D(i, j) + lambda_E(i, j)) / 4.0;
                        double a_ii = C * h * h / tau + GAMMA * sum_L;
                        uh_new[m] = Fh[m] / a_ii; // D^{-1} * Fh
                        if (is_internal(i + 1, j)) { // 右側
                            uh_new[m] += GAMMA * uh_old[indices[i + 1][j]] * L_B / a_ii;
                        }
                        if (is_internal(i - 1, j)) { // 左側
                            uh_new[m] += GAMMA * uh_old[indices[i - 1][j]] * L_C / a_ii;
                        }
                        if (is_internal(i, j + 1)) { // 上側
                            uh_new[m] += GAMMA * uh_old[indices[i][j + 1]] * L_D / a_ii;
                        }
                        if (is_internal(i, j - 1)) { // 下側
                            uh_new[m] += GAMMA * uh_old[indices[i][j - 1]] * L_E / a_ii;
                        }
                    }
                }
            }

            // 未知ベクトルを更新する
            double uh_norm = 0.0;
            double res = 0.0;
            for (int m = 0; m < n; ++m) {
                res += (uh_old[m] - uh_new[m]) * (uh_old[m] - uh_new[m]);
                uh_norm += uh_old[m] * uh_old[m];
                uh_old[m] = uh_new[m];
            }
            res /= uh_norm;
            // 更新分の二乗和を見て反復計算を続けるか判断する
            if (res < TOLERANCE || l == MAXITER - 1) {
                // 標準エラー出力にログを書き出す．
                // このようにデバッグ用のログは標準エラー出力に出しておけば，
                // 標準出力にでた数値結果の方を汚さずにファイルに保存できる．
//                fprintf(stderr, "[k = %d][%d] res = %e\n",k, l, res);
                break;
            }
        }
        uh[k + 1] = uh_new;

    }

    // 出力する
    double gosa = 0.0; // error
    double max_gosa = 0.0;
    int vec_len = 0; // vector length
    for (int k = 1; k <= r; k++) {
        for (int i = 0; i <= N; i += CULL) {
            for (int j = 0; j <= N; j += CULL) {
                if (is_internal(i, j)) {
                    int m = indices[i][j];
//                    if(k==r){
//                        printf("%f %f %f\n", X0 + i * h, Y0 + j * h, uh[k][m]);
//                    }
                    // debug
//                    if (fabs(X0 + i * h) < h / 8 && fabs(Y0 + j * h) < h / 8) {
//                        fprintf(stderr, "U(0, 0, %.2f) = %.2e\t u = %.2e\n", k2t(k), uh[k][m], Dirichlet_data(i2x(i), j2y(j),k2t(k)));
//                    }
                    // error
                    gosa += pow(uh[k][m] - Dirichlet_data(i2x(i), j2y(j), k2t(k)), 2);
                    vec_len++;
                }

            }
        }
        if(gosa>max_gosa){
            max_gosa= gosa;

        }
        gosa = 0.0;
    }
//    double gosa_norm = sqrt(gosa) / vec_len;
    double gosa_norm = sqrt(max_gosa) * sqrt(h);
    fprintf(stderr, "N = %d, r = %d, error = %.2e\n", N, r, gosa_norm);
    fprintf(fr, "%d %d %e %e %e\n", N,r, h,tau, gosa_norm);

    // // // // // // // // // // * 重 要 * // // // // // // // // // //
    // // // // // malloc で確保したメモリは必ず解放する！ // // // // //
    for (int k = 0; k <= r; k++) {
        free(uh[k]);
    }
    free(uh);
    free(Fh);
    fclose(fd);
    for (int i = 0; i <= N; i++) {
        free(indices[i]);
    }
    free(indices);
}

// // // // // // // // // * キ ケ ン * // // // // // // // // // //
// // // // // malloc 関数の使用には十分注意すること！ // // // // //
// 事前に配列サイズが分からないときは，メモリ確保命令 malloc を使い
// 浮動小数点数型の配列を動的確保する．
// 使い終わったら必ず free() すること．
double *allocate_real_vector(size_t n) {
    double *vector = malloc(n * sizeof(double));
    if (vector == NULL) {
        fprintf(stderr, "Failed to allocate vector. ABORT.\n");
        exit(EXIT_FAILURE);
    }
    return vector;
}

