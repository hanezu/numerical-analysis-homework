// neumann condition on the upper hand circle boundary
// k = 2, phi_N = 4*(x^2 + y^2)
// by: Huachun Zhu


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include "poisson-neumann-region.h"

// 矩形領域 [X0, X0 + L] \times [Y0, Y0 + L] に格子点を作る
// 格子の基準点
#define X0 (-1.0)
#define Y0 (-1.0)
// 格子の縦および横の長さ
#define L (2.0)
// L あたりの分割数
int N = 50;
// 分割幅
#define h (L / N)

// ヤコビ法の反復終了の閾値
#define TOLERANCE (1.0e-13)
// ヤコビ法の最大反復回数
#define MAXITER (20000)
// 出力時に間引くデータ間隔
#define CULL (4)

// 0 if regular, 1 if near dirichlet, 2 if near neumann, -1 if external
#define REGULAR (0)
#define DIRICHLET (1)
#define NEUMANN (2)
#define EXTERNAL (-1)

FILE *fd;
FILE *fr;


// change indexes to coordinates
double i2x(int i) {
    return X0 + i * h;
}

double j2y(int j) {
    return Y0 + j * h;
}

// put fundamental vector stuffs here
typedef struct Vector {
    double x;
    double y;
} Vector;

typedef struct Point {
    int i;
    int j;
} Point;

typedef struct X_with_lambda {
    Vector *X;
    Point *Y;
    Point *Z;
    double lambda_Y;
    double lambda_Z;
} X_with_lambda;


Vector *vector_ij(int i, int j) {
    // a simple ij constructor
    Vector *v = malloc(sizeof(Vector));
    v->x = i2x(i);
    v->y = j2y(j);
    return v;
}

Point *point_ij(int i, int j) {
    // a simple ij constructor
    Point *p = malloc(sizeof(Point));
    p->i = i;
    p->j = j;
    return p;
}

Vector *vector_xy(double x, double y) {
    // a simple ij constructor
    Vector *v = malloc(sizeof(Vector));
    v->x = x;
    v->y = y;
    return v;
}

//
//Vector *add(Vector *v1, Vector *v2){
//    Vector *v = malloc(sizeof(Vector));
//    v->x = v1->x + v2->x;
//    v->y = v1->y + v2->y;
//    return v;
//}
//
//Vector *subtract(Vector *v1, Vector *v2){
//    Vector *v = malloc(sizeof(Vector));
//    v->x = v1->x - v2->x;
//    v->y = v1->y - v2->y;
//    return v;
//}
//
void magnify(Vector *v, double a) {
    v->x *= a;
    v->y *= a;
}

double length_vec(Vector *vector) {
    return sqrt(vector->x * vector->x + vector->y * vector->y);
}

void normalize(Vector *vector) {
    double vector_len = length_vec(vector);
    vector->x /= vector_len;
    vector->y /= vector_len;
}


// 内部節点に連番をふるための配列
//int indices[N + 1][N + 1];
int **indices;

// ポアソン方程式の右辺
double f(double x, double y) {
    return 0.0;
}

// 境界条件
double Dirichlet_data(double x, double y) {
    return x * x - y * y;
}

double Neumann_k(double x, double y) {
    return 2.0;
}

double Neumann_phi(double x, double y) {
    return 4 * (x * x - y * y);
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

// check if a boundary point near neumann. false means near dirichlet
bool near_neumann(int i, int j) {
    // you can change the value for different neumann boundary
    // if set it to -1, there will be no dirichlet boundary
    if (Y0 + j * h >= 0) {
        //debug for on the boundary
//        double x = i2x(i);
//        double y = j2y(j);
//        if(fmin(fabs(y-y_top(x)),fabs(y-y_bottom(x)))+fmin(fabs(x-x_right(y)),fabs(x-x_left(y)))==0){
//            fprintf(stderr, "(%d, %d) is on the neumann boundary!", i,j);
//        }
        return true;
    } else {
        return false;
    }
}

double avg_lambda(int i, int j);

//　return the condition of point
int point_condition(int i, int j) {
    if (is_internal(i, j)) {
        double C = avg_lambda(i, j);
        if (C == 1) {
            return REGULAR;
        } else {
            if (near_neumann(i, j)) {
                return NEUMANN;
            } else {
                return DIRICHLET;
            }
        }
    } else {
        return EXTERNAL;
    }

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


//find the point of vertical
Vector *vertical(int i, int j) {
    //calc vector (x0,y0), and normalize it
//    double x0 = i2x(i);
//    double y0 = j2y(j);
//    Vector *from_0_0 = malloc(sizeof(Vector));
//    from_0_0 = &(Vector){
//        .x = x0, .y = y0
//    };
    Vector *from_0_0 = vector_ij(i, j);
    normalize(from_0_0);
    return from_0_0;

}

//find the point cross at back
X_with_lambda *cross_backward(int i, int j) {
    // if vector has |y|>|x|, it will cross back at the j-1 line
    Vector *direction_vector;
    Vector *point;
    Point *Y;
    Point *Z;
    double x;
    double y;
    double magnify_ratio;
    double lambda_Y;
    double lambda_Z;

    point = vector_ij(i, j);
    direction_vector = vertical(i, j);
    x = direction_vector->x;
    y = direction_vector->y;
//    fprintf(fd,"current point: (%3d, %3d)\t",i,j);
//    assert(y>=0);

    if (fabs(y) >= fabs(x)) {
        //|k|>=1, cross on the lower line
        // magnify by h/|y|
        magnify_ratio = 1 / fabs(y);
        magnify(direction_vector, magnify_ratio);
        if (y >= 0) {
            Y = point_ij(i, j - 1);
            if (x > 0) {
                Z = point_ij(i - 1, j - 1);
            } else {
                Z = point_ij(i + 1, j - 1);
            }
        } else {
            Y = point_ij(i, j + 1);
            if (x > 0) {
                Z = point_ij(i - 1, j + 1);
            } else {
                Z = point_ij(i + 1, j + 1);
            }
        }

        lambda_Y = fabs(direction_vector->x);
        lambda_Z = 1 - lambda_Y;

    } else {
        // |k|<=1, cross on the vertical line
        // magnify by h/|x|
        // the relation between Y and Z? don't know
        // but it doesn't matter because of the symmetrical property between X and Y
        magnify_ratio = 1 / fabs(x);
        magnify(direction_vector, magnify_ratio);
        if (x >= 0) {
            Y = point_ij(i - 1, j);
            if (y >= 0) {
                Z = point_ij(i - 1, j - 1);
            } else {
                Z = point_ij(i - 1, j + 1);
            }
        } else {
            Y = point_ij(i + 1, j);
            if (y >= 0) {
                Z = point_ij(i + 1, j - 1);
            } else {
                Z = point_ij(i + 1, j + 1);
            }
        }

        lambda_Y = fabs(direction_vector->x);
        lambda_Z = 1 - lambda_Y;
    }
    magnify(direction_vector, h);

//    fprintf(fd,"Y: (%3d, %3d)\t",Y->i,Y->j);
//    fprintf(fd,"Z: (%3d, %3d)\n",Z->i,Z->j);
    x = direction_vector->x;
    y = direction_vector->y;
    point->x -= x;
    point->y -= y;

    X_with_lambda *x_lambda = malloc(sizeof(X_with_lambda));
    x_lambda->X = point;
    x_lambda->Y = Y;
    x_lambda->Z = Z;
    x_lambda->lambda_Y = lambda_Y;
    x_lambda->lambda_Z = lambda_Z;
    return x_lambda;


}

// この関数は下で定義する．
double *allocate_real_vector(size_t);

void test(void) {
    int i = 108;
    int j = 108;
    Vector *v = vector_ij(i, j);
    fprintf(stderr, "%lf, %lf\n", v->x, v->y);
    Vector *dv = vertical(i, j);
    fprintf(stderr, "%lf, %lf\n", dv->x, dv->y);
    X_with_lambda *bv = cross_backward(i, j);
    fprintf(stderr, "%lf, %lf\n", bv->X->x, bv->X->y);

}

void run_main(int split_times);

int main(int argc, char *argv[]) {
    bool is_test = false;
    if (is_test) {
        test();
    } else {
//        int N_array[] = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190};
//        int N_array[] = {100};
        int N_array[] = {10,20,30,40,50,60,70,80,90};

        int len_N_array = sizeof(N_array)/ sizeof(N_array[0]);
        fr = fopen("error_fit.csv", "w+");
        fprintf(fr,"iter_times N h error\n");
        for(int i = 0; i < len_N_array; i++){
            run_main(N_array[i]);
        }
        fclose(fr);
    }

    return 0;
}

void run_main(int split_times) {
    N = split_times;
//    fd = fopen("debug.txt", "w");

    indices = malloc((N + 1) * sizeof(*indices)); //init indices row

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
    double *uh = allocate_real_vector(n);
    double *uh_new = allocate_real_vector(n);
    double *Fh = allocate_real_vector(n);

    // 右辺ベクトルを初期化する
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            if (is_internal(i, j)) {
                int m = indices[i][j];
                double x = X0 + i * h;
                double y = Y0 + j * h;
                int cond = point_condition(i, j);
                if (cond <= DIRICHLET) {//dirichlet or regular
                    // 正則な節点の場合，C = 1 となることに注意
                    double C = (lambda_B(i, j) + lambda_C(i, j)) * (lambda_D(i, j) + lambda_E(i, j)) / 4.0;
                    Fh[m] = C * f(x, y) * h * h;
                    // 境界条件から既知の数を右辺ベクトルに足し上げる．
                    if (!is_internal(i + 1, j)) { // 右側
                        Fh[m] += Dirichlet_data(x_right(y), y) / (lambda_B(i, j) * h * h);
                    }
                    if (!is_internal(i - 1, j)) { // 左側
                        Fh[m] += Dirichlet_data(x_left(y), y) / (lambda_C(i, j) * h * h);
                    }
                    if (!is_internal(i, j + 1)) { // 上側
                        Fh[m] += Dirichlet_data(x, y_top(x)) / (lambda_D(i, j) * h * h);
                    }
                    if (!is_internal(i, j - 1)) { // 下側
                        Fh[m] += Dirichlet_data(x, y_bottom(x)) / (lambda_E(i, j) * h * h);
                    }
//                    if(C==1){
//                        fprintf(fd,"Regular: Fh[%3d] = %.3lf\n", m,Fh[m]);
//                    }else{
//                        fprintf(fd,"Dirichlet: Fh[%3d] = %.3lf\n", m,Fh[m]);
//                    }
                } else if (cond == NEUMANN) {
                    Vector *Q = vector_ij(i, j);
                    normalize(Q);
                    Fh[m] = Neumann_phi(Q->x, Q->y);
//                    fprintf(fd,"Neumann: Fh[%3d] = %.3lf\n", m,Fh[m]);
                }
            }
        }
    }

    // ヤコビ法による反復計算で方程式 Ah * uh = Fh を解く
    // 反復の初期値は D^{-1} * Fh に近いものでよい
    for (int m = 0; m < n; ++m) {
        uh[m] = Fh[m] * h * h / 4.0;
    }

    for (int k = 0; k < MAXITER; ++k) {
        // 疎行列 H = -D^{-1} * R は非零要素を非対角に高々 4 つもつ
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                if (is_internal(i, j)) {
                    int cond = point_condition(i, j);
                    int m = indices[i][j];
                    if (cond <= DIRICHLET) { // dirichlet or regular
                        double L_B = 1.0 / lambda_B(i, j);
                        double L_C = 1.0 / lambda_C(i, j);
                        double L_D = 1.0 / lambda_D(i, j);
                        double L_E = 1.0 / lambda_E(i, j);
                        double a_ii = L_B + L_C + L_D + L_E;
                        uh_new[m] = Fh[m] * h * h / a_ii; // D^{-1} * Fh
                        if (is_internal(i + 1, j)) { // 右側
                            uh_new[m] += uh[indices[i + 1][j]] * L_B / a_ii;
                        }
                        if (is_internal(i - 1, j)) { // 左側
                            uh_new[m] += uh[indices[i - 1][j]] * L_C / a_ii;
                        }
                        if (is_internal(i, j + 1)) { // 上側
                            uh_new[m] += uh[indices[i][j + 1]] * L_D / a_ii;
                        }
                        if (is_internal(i, j - 1)) { // 下側
                            uh_new[m] += uh[indices[i][j - 1]] * L_E / a_ii;
                        }
                    } else if (cond == NEUMANN) {
                        Vector *A = vector_ij(i, j);
                        X_with_lambda *x_with_lambda = cross_backward(i, j);
                        Vector *X = x_with_lambda->X;
                        Point *Y = x_with_lambda->Y;
                        Point *Z = x_with_lambda->Z;
                        Vector *Q = vector_ij(i, j);
                        normalize(Q);
                        double distance_AX = sqrt(pow(A->x - X->x, 2.0) + pow(A->y - X->y, 2.0));
                        double a_ii = (1 + distance_AX * Neumann_k(Q->x, Q->y)) / distance_AX;
                        uh_new[m] = Fh[m] / a_ii; // D^{-1} * Fh

                        // D^{-1} * {A[Y]*u[Y]+A[Z]*u[Z]
                        uh_new[m] += uh[indices[Y->i][Y->j]] * x_with_lambda->lambda_Z / distance_AX / a_ii;
                        uh_new[m] += uh[indices[Z->i][Z->j]] * x_with_lambda->lambda_Y / distance_AX / a_ii;
//                        fprintf(fd,"at A(%d, %d) = (%.3lf, %.3lf), X(%.3lf, %.3lf), Y(%d, %d), Z(%d, %d), Q(%.3lf, %.3lf), ", i,j,A->x,A->y,X->x,X->y,Y->i,Y->j,Z->i,Z->j,Q->x,Q->y);
//                        fprintf(fd,"|A-X| = %.3lf, lambda_Y = %.3lf, lambda_Z = %.3lf, uh = %.3lf, uh_new = %.3lf\n",distance_AX,x_with_lambda->lambda_Y, x_with_lambda->lambda_Z,uh[m],uh_new[m]);

                    }
                }
            }
        }

        // 未知ベクトルを更新する
        double uh_norm = 0.0;
        double res = 0.0;
        for (int m = 0; m < n; ++m) {
            res += (uh[m] - uh_new[m]) * (uh[m] - uh_new[m]);
            uh_norm += uh[m] * uh[m];
            uh[m] = uh_new[m];
        }
        res /= uh_norm;
        // 更新分の二乗和を見て反復計算を続けるか判断する
        if (res < TOLERANCE || k == MAXITER - 1) {
            // 標準エラー出力にログを書き出す．
            // このようにデバッグ用のログは標準エラー出力に出しておけば，
            // 標準出力にでた数値結果の方を汚さずにファイルに保存できる．
            fprintf(stderr, "[%d] res = %e\n", k, res);
            fprintf(fr,"%d ",k);
            break;
        }
    }

    // 出力する
    double gosa = 0.0; // error
    int vec_len = 0; // vector length
    for (int i = 0; i <= N; i += CULL) {
        for (int j = 0; j <= N; j += CULL) {
            if (is_internal(i, j)) {
                int m = indices[i][j];
                printf("%f %f %f\n", X0 + i * h, Y0 + j * h, uh[m]);
                // debug
                if (fabs(X0 + i * h) < h / 8 && fabs(Y0 + j * h) < h / 8) {
                    fprintf(stderr, "u(0, 0) = %e\n", uh[m]);
                }
                // error
                gosa += pow(uh[m] - Dirichlet_data(i2x(i), j2y(j)), 2);
                vec_len++;
            }

        }
    }
//    double gosa_norm = sqrt(gosa) / vec_len;
    double gosa_norm = sqrt(gosa) * sqrt(h);
    fprintf(stderr, "N = %d\terror = %e\n", N, gosa_norm);
    fprintf(fr, "%d %e %e\n",N,h, gosa_norm);

    // // // // // // // // // // * 重 要 * // // // // // // // // // //
    // // // // // malloc で確保したメモリは必ず解放する！ // // // // //
    free(uh);
    free(uh_new);
    free(Fh);
//    fclose(fd);
    for (int i = 0;i<=N;i++){
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

