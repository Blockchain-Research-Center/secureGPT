#include <chrono>
#include <cstdio>
#include <future>
#include <iostream>
#include <memory>
#include <random>
#include <seal/ciphertext.h>
#include <seal/evaluator.h>
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
    int N = 1 << 12, depth = 19;
    double init_scale = pow(2.0, 40);
    shared_ptr<CKKSManager> ckks = make_shared<CKKSManager>(N, depth, init_scale, 1);
    SgnEvaluator sgn_eval(0, ckks->encoder, 0.5);

    vector<double> x{ 4.2, 3.9, 1.8, -2.5, -8.6, 6 };
    // vector<double> pt = init_vec(2000);
    Ciphertext ct, b0, b1, b2;
    Plaintext p0, p1, p2, delta;
    vector<double> dest;
    ckks->encode_and_encrypt(x, ct);

    ckks->encoder->encode(init_vec_with_value(ckks->real_number_capacity(), -4.0), ct.parms_id(), ct.scale(), p0);
    ckks->encoder->encode(init_vec_with_value(ckks->real_number_capacity(), -1.95), ct.parms_id(), ct.scale(), p1);
    ckks->encoder->encode(init_vec_with_value(ckks->real_number_capacity(), 3.0), ct.parms_id(), ct.scale(), p2);
    ckks->encoder->encode(
        init_vec_with_value(ckks->real_number_capacity(), 1.0 / 16), ct.parms_id(), ct.scale(), delta);

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
    ckks->encoder->encode(
        init_vec_with_value(ckks->real_number_capacity(), 0.5), b0.parms_id(), b0.scale(), zero_point_five);
    Ciphertext a0, a1, a2, a3;
    // Ciphertext neg_b0;
    // ckks->evaluator->negate(b0,neg_b0);
    // ckks->evaluator->add_plain(neg_b0, zero_point_five, a0); // a0 = -b0 + 0.5
    ckks->evaluator->sub(b0, b1, a1);                    // a1 = b0 - b1
    ckks->evaluator->sub(b1, b2, a2);                    // a2 = b1 - b2
    ckks->evaluator->add_plain(b2, zero_point_five, a3); // a3 = b2 + 0.5

    time_end = high_resolution_clock::now();
    cout << ckks->real_number_capacity() << " times GELU" << endl;
    cout << "Computation cost:  " << duration_cast<milliseconds>(time_end - time_start).count() / 2 << " ms"
         << endl; // TODO: HEXL Acc
    cout << "Communication cost:  " << a0.save_size(compr_mode_type::zstd) << " bytes" << endl;

    //============================================
    double A[] = { -0.5054031199708174, -0.4222658115198386, -0.1180761295118195, -0.0110341340306157 };
    double B[] = { 0.5, 0.0085263215410380, 0.3603292692789629, 0, -0.037688200365904, 0, 0.0018067462606141 };

    Ciphertext px;
    ckks->evaluator->sub(ct, ct, px);

    // cout << a1.scale() << " " << a2.scale() << " " << a3.scale() << endl;

    uint64_t p = ckks->get_modulus(ct, 1);
    uint64_t q = ckks->get_modulus(ct, 2);
    uint64_t r = ckks->get_modulus(ct, 3);
    uint64_t s = ckks->get_modulus(ct, 4);
    uint64_t t = ckks->get_modulus(ct, 5);

    uint64_t P = ckks->get_modulus(a1, 1);

    double D = init_scale / a1.scale() * P;

    Ciphertext x2, x3, x6;

    ckks->evaluator->square(ct, x2);
    ckks->evaluator->relinearize_inplace(x2, ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(x2); // L-1

    auto L = ct.parms_id(); // L

    Ciphertext ct_next;
    ckks->evaluator->mod_reduce_to_next(ct, ct_next);
    // ckks->evaluator->mod_switch_to_next_inplace(ct); // L-1

    ckks->evaluator->multiply(x2, ct_next, x3);
    ckks->evaluator->relinearize_inplace(x3, ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(x3); // L-2

    ckks->evaluator->square(x3, x6);
    ckks->evaluator->relinearize_inplace(x6, ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(x6); // L-3

    Plaintext A0, A1, A2, A3;
    double a2_scale = D / x2.scale() * q;
    ckks->encoder->encode(A[0], x3.parms_id(), D, A0);
    ckks->encoder->encode(A[1], ct_next.parms_id(), D / ct_next.scale() * q, A1);
    ckks->encoder->encode(A[2], x2.parms_id(), a2_scale, A2);
    ckks->encoder->encode(A[3], ct.parms_id(), a2_scale / ct.scale() * p, A3);

    Ciphertext a1x;
    ckks->evaluator->multiply_plain(ct_next, A1, a1x);
    ckks->evaluator->rescale_to_next_inplace(a1x);

    ckks->evaluator->add_plain_inplace(a1x, A0);

    Ciphertext a3x;
    ckks->evaluator->multiply_plain(ct, A3, a3x);
    ckks->evaluator->rescale_to_next_inplace(a3x);

    ckks->evaluator->add_plain_inplace(a3x, A2);

    ckks->evaluator->multiply_inplace(a3x, x2);
    ckks->evaluator->relinearize_inplace(a3x, ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(a3x);

    ckks->evaluator->add_inplace(a3x, a1x);

    Plaintext B0, B1, B2, B3, B4, B5, B6;

    ckks->encoder->encode(B[0], x3.parms_id(), D, B0);
    ckks->encoder->encode(B[1], ct_next.parms_id(), D / ct_next.scale() * q, B1);
    ckks->encoder->encode(B[2], x2.parms_id(), D / x2.scale() * q, B2);

    double a3_scale = D / x3.scale() * r;
    ckks->encoder->encode(B[3], x3.parms_id(), a3_scale, B3);
    ckks->encoder->encode(B[4], ct_next.parms_id(), a3_scale / ct_next.scale() * q, B4);
    ckks->encoder->encode(B[5], x2.parms_id(), a3_scale / x2.scale() * q, B5);

    ckks->encoder->encode(B[6], x6.parms_id(), D / x6.scale() * s, B6);

    Ciphertext b1x;
    ckks->evaluator->multiply_plain(ct_next, B1, b1x);
    ckks->evaluator->rescale_to_next_inplace(b1x);

    ckks->evaluator->add_plain_inplace(b1x, B0);

    Ciphertext b2x2;
    ckks->evaluator->multiply_plain(x2, B2, b2x2);
    ckks->evaluator->rescale_to_next_inplace(b2x2);

    ckks->evaluator->add_inplace(b2x2, b1x);

    Ciphertext b4x;
    ckks->evaluator->multiply_plain(ct_next, B4, b4x);
    ckks->evaluator->rescale_to_next_inplace(b4x);

    ckks->evaluator->add_plain_inplace(b4x, B3);

    Ciphertext b5x2;
    ckks->evaluator->multiply_plain(x2, B5, b5x2);
    ckks->evaluator->rescale_to_next_inplace(b5x2);

    ckks->evaluator->add_inplace(b5x2, b4x);

    ckks->evaluator->multiply_inplace(b5x2, x3);
    ckks->evaluator->relinearize_inplace(b5x2, ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(b5x2);

    Ciphertext b6x6;
    ckks->evaluator->multiply_plain(x6, B6, b6x6);
    ckks->evaluator->rescale_to_next_inplace(b6x6);

    ckks->evaluator->mod_switch_to_inplace(b2x2, b6x6.parms_id());
    ckks->evaluator->mod_switch_to_inplace(b5x2, b6x6.parms_id());

    ckks->evaluator->add_inplace(b5x2, b2x2);
    ckks->evaluator->add_inplace(b6x6, b5x2);

    Ciphertext s1, s2, s3;
    Ciphertext new_ct;
    Plaintext new_pt;

    ckks->encoder->encode(x, a3x.scale(), new_pt);
    ckks->encryptor->encrypt(new_pt, new_ct);

    ckks->evaluator->mod_switch_to_inplace(a3x, a1.parms_id());
    ckks->evaluator->multiply(a3x, a1, s1);
    ckks->evaluator->relinearize_inplace(s1, ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(s1);

    ckks->evaluator->mod_switch_to_inplace(b6x6, a2.parms_id());
    ckks->evaluator->multiply(b6x6, a2, s2);
    ckks->evaluator->relinearize_inplace(s2, ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(s2);

    ckks->evaluator->mod_switch_to_inplace(new_ct, a3.parms_id());
    ckks->evaluator->multiply(new_ct, a3, s3);
    ckks->evaluator->relinearize_inplace(s3, ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(s3);

    ckks->evaluator->add_inplace(s2, s1);
    ckks->evaluator->add_inplace(s3, s2);

    //============================================

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

    ckks->decrypt_and_decode(s3, dest);
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
