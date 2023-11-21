#include <chrono>
#include <cstdio>
#include <future>
#include <iostream>
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
    int N = 1 << 12, depth = 17;
    double init_scale = pow(2.0, 40);
    shared_ptr<CKKSManager> ckks = make_shared<CKKSManager>(N, depth, init_scale, 1);
    SgnEvaluator sgn_eval(0, ckks->encoder, 0.5);

    vector<double> x{4.2, 3.9, 1.8, -2.5, -8.6};
    // vector<double> pt = init_vec(2000);
    Ciphertext ct, b0, b1, b2;
    Plaintext p0, p1, p2, delta;
    vector<double> dest;
    ckks->encode_and_encrypt(x, ct);

    ckks->encoder->encode(init_vec_with_value(ckks->real_number_capacity(), -4.0), ct.parms_id(), ct.scale(), p0);
    ckks->encoder->encode(init_vec_with_value(ckks->real_number_capacity(), -1.95), ct.parms_id(), ct.scale(), p1);
    ckks->encoder->encode(init_vec_with_value(ckks->real_number_capacity(), 3.0), ct.parms_id(), ct.scale(), p2);
    ckks->encoder->encode(init_vec_with_value(ckks->real_number_capacity(), 1.0/16), ct.parms_id(), ct.scale(), delta);

    ckks->evaluator->sub_plain(ct, p0, b0);
    ckks->evaluator->multiply_plain_inplace(b0, delta);
    ckks->evaluator->rescale_to_next_inplace(b0);
    ckks->evaluator->sub_plain(ct, p1, b1);
    ckks->evaluator->multiply_plain_inplace(b1, delta);
    ckks->evaluator->rescale_to_next_inplace(b1);
    ckks->evaluator->sub_plain(ct, p2, b2);
    ckks->evaluator->multiply_plain_inplace(b2, delta);
    ckks->evaluator->rescale_to_next_inplace(b2);

    time_start = high_resolution_clock::now();
    sgn_eval.sgn(2, 2, b0, b0, ckks);
    sgn_eval.sgn(2, 2, b1, b1, ckks);
    sgn_eval.sgn(2, 2, b2, b2, ckks);

    Plaintext zero_point_five;
    ckks->encoder->encode(init_vec_with_value(ckks->real_number_capacity(), 0.5), b0.parms_id(), b0.scale(), zero_point_five);
    Ciphertext a0, a1, a2, a3;
    // Ciphertext neg_b0;
    // ckks->evaluator->negate(b0,neg_b0);
    // ckks->evaluator->add_plain(neg_b0, zero_point_five, a0); // a0 = -b0 + 0.5
    ckks->evaluator->sub(b0, b1, a1); // a1 = b0 - b1
    ckks->evaluator->sub(b1, b2, a2); // a2 = b1 - b2
    ckks->evaluator->add_plain(b0, zero_point_five, a3); // a3 = b2 + 0.5


    time_end = high_resolution_clock::now();
    cout  << ckks->real_number_capacity() << " times GELU" << endl;
    cout << "Computation cost:  " << duration_cast<milliseconds>(time_end - time_start).count() / 2 << " ms"
         << endl; // TODO: HEXL Acc
    cout << "Communication cost:  " << a0.save_size(compr_mode_type::zstd) << " bytes" << endl;

    
    // uint64_t t = ckks->get_modulus(ct, 5);

    // ckks->evaluator->mod_reduce_to_inplace(ct, score.parms_id());
    // ct.scale() = t;

    // ckks->evaluator->multiply(score, ct, relu_x);
    // ckks->evaluator->relinearize_inplace(relu_x, ckks->relin_keys);
    // ckks->evaluator->rescale_to_next_inplace(relu_x);

    // time_end = high_resolution_clock::now();
    // cout  << N << " times Relu computation" << endl;
    // cout << "Computation cost:  " << duration_cast<milliseconds>(time_end - time_start).count() / 2 << " ms"
    //      << endl; // TODO: HEXL Acc
    // cout << "Communication cost:  " << relu_x.save_size(compr_mode_type::zstd) << " bytes" << endl;
    ckks->decrypt_and_decode(a1, dest);
    // cout << "input: ";
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
