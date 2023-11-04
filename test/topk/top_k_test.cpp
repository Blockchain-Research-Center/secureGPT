#include <chrono>
#include <future>
#include <memory>
#include <random>
#include <seal/plaintext.h>
#include <seal/seal.h>
#include <seal/valcheck.h>
#include <thread>
#include <vector>
#include "sgn.h"

using namespace std;
using namespace seal;
using namespace std::chrono;

vector<double> init_vec(int N);
vector<double> init_vec_with_value(int N, double init_value);

int main()
{
    chrono::high_resolution_clock::time_point time_start, time_end;
    int N = 4096;
    double init_scale = pow(2.0, 50);
    int dg1 = 6;      // Argmax 1st polynomial (dg)
    int df1 = 1;      // Argmax 1st polynomial (df), set to 0 for binary
    int dg2 = 2;      // Argmax 2nd polynomial (dg), set to 0 for binary
    int df2 = 2;      // Argmax 2nd polynomial (df)
    int depth_nn = 5; // 1 for linear, 2 for learnable act, 1 before linear for rot
    int depth_reduce = 1;
    int depth_argmax = 4 * (dg1 + df1 + dg2 + df2); // potentially add one to mask out between sgns
    int depth_mask_winner = 1;
    // int depth = depth_nn + depth_reduce + depth_argmax + depth_mask_winner;
    int depth = 16 + 1;
    cout << "depth: " << depth << endl;
    shared_ptr<CKKSManager> ckks = make_shared<CKKSManager>(N, depth, init_scale, 1);
    SgnEvaluator sgn_eval(0, ckks->encoder, 0.5);

    vector<double> x{ 0.6, 0.4, -0.8, 0.2, -0.9, 0.5 };
    // vector<double> pt = init_vec(2000);
    Ciphertext ct;
    Ciphertext score;
    Ciphertext relu_x;
    vector<double> dest;
    ckks->encode_and_encrypt(x, ct);
    time_start = high_resolution_clock::now();
    sgn_eval.sgn(dg2, df2, ct, score, ckks);

    Plaintext plain_one;
    ckks->encoder->encode(init_vec_with_value(x.size(), 0.5), score.parms_id(), score.scale(), plain_one);

    ckks->evaluator->add_plain_inplace(score, plain_one);

    uint64_t t = ckks->get_modulus(ct, 5);

    ckks->evaluator->mod_reduce_to_inplace(ct, score.parms_id());
    ct.scale() = t;

    ckks->evaluator->multiply(score, ct, relu_x);
    ckks->evaluator->relinearize_inplace(relu_x, ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(relu_x);

    time_end = high_resolution_clock::now();
    cout << "SGN took " << duration_cast<milliseconds>(time_end - time_start).count() / 2 << " ms"
         << endl; // TODO: HEXL Acc
    cout << "Comm  " << score.save_size(compr_mode_type::zstd) << " " << endl;
    ckks->decrypt_and_decode(relu_x, dest);
    cout << "result: ";
    for (int i = 0; i < x.size(); i++) {
        cout << dest[i] << " ";
    }
    cout << endl;
}

vector<double> init_vec(int N)
{
    std::vector<double> v(N);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    for (int i = 0; i < N; ++i) {
        v[i] = dis(gen);
    }

    return v;
}

vector<double> init_vec_with_value(int N, double init_value)
{
    std::vector<double> v(N);

    for (int i = 0; i < N; ++i) {
        v[i] = init_value;
    }

    return v;
}
