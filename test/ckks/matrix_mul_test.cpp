#include <chrono>
#include <cstdlib>
#include <memory>
#include <random>
#include <seal/ciphertext.h>
#include <seal/evaluator.h>
#include <seal/plaintext.h>
#include <seal/seal.h>
#include <seal/util/polyarithsmallmod.h>
#include <seal/valcheck.h>
#include <thread>
#include <vector>
#include "ckks_manager.h"

using namespace std;
using namespace seal;
using namespace std::chrono;
using namespace seal::util;

shared_ptr<CKKSManager> ckks;

vector<double> gen_random_vector(size_t n, mt19937 &gen);
vector<Ciphertext> expand_ciphertext(
    const Ciphertext &encrypted, uint32_t m, GaloisKeys &galkey, vector<int> &galois_elts);
void multiply_power_of_X(const Ciphertext &encrypted, Ciphertext &destination, uint32_t index);

int main()
{
    chrono::high_resolution_clock::time_point time_start, time_end;
    mt19937 gen(random_device{}());
    int N = 1 << 12;
    double init_scale = pow(2.0, 40);
    int depth = 1;
    ckks = make_shared<CKKSManager>(N, depth, init_scale, 1);

    // {
    //     Plaintext pt;
    //     auto scale = pow(2.0, 50);
    //     vector<double> v = { 1.0, 2.0 };

    //     ckks->encoder->encode(v, scale, pt);
    //     // std::cout << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;

    //     Ciphertext ct;
    //     ckks->encryptor->encrypt(pt, ct);
    //     ckks->evaluator->rotate_vector_inplace(ct, -1, ckks->galois_keys);
    //     ckks->decryptor->decrypt(ct, pt);
    //     ckks->encoder->decode(pt, v);
    //     for (auto &e : v) {
    //         std::cout << e << " ";
    //     }
    //     // std::cout << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
    // }

    // exit(0);

    vector<Plaintext> a_pts;
    a_pts.reserve(768);
    for (int i = 0; i < 768; i++) {
        Plaintext pt;
        ckks->encoder->encode(gen_random_vector(ckks->real_number_capacity(), gen), init_scale, pt);
        a_pts.emplace_back(pt);
    }

    vector<Ciphertext> b_compressed_cts;
    for (int i = 0; i < 768 * 768 / ckks->real_number_capacity(); i++) {
        Plaintext pt;
        ckks->encoder->encode(gen_random_vector(ckks->real_number_capacity(), gen), init_scale, pt);
        Ciphertext ct;
        ckks->encryptor->encrypt(pt, ct);
        b_compressed_cts.push_back(ct);
    }

    time_start = high_resolution_clock::now();
    vector<Ciphertext> b_expanded_cts;
    for (Ciphertext &ct : b_compressed_cts) {
        vector<Ciphertext> temp_cts = expand_ciphertext(ct, ckks->poly_modulus_degree(), ckks->galois_keys, ckks->rots);
        cout << "expanding..." << endl;
        b_expanded_cts.insert(
            b_expanded_cts.end(), make_move_iterator(temp_cts.begin()), make_move_iterator(temp_cts.end()));
    }
    time_end = high_resolution_clock::now();
    cout << "expanding time: " << duration_cast<seconds>(time_end - time_start).count() << "seconds" << endl;

    vector<Ciphertext> res_cts;
    Plaintext pt(ckks->poly_modulus_degree(), 0);
    time_start = high_resolution_clock::now();
    for (int i = 0; i < 768; i++) {
        Ciphertext res_col_ct;
        ckks->encryptor->encrypt(pt, res_col_ct);
        for (int j = 0; j < 768; j++) {
            Ciphertext temp;
            ckks->evaluator->multiply_plain(b_expanded_cts[i * 768 + j], a_pts[j], temp);
            ckks->evaluator->add(res_col_ct, temp, res_col_ct);
        }
        res_cts.push_back(res_col_ct);
    }
    time_end = high_resolution_clock::now();
    cout << "calculating res time: " << duration_cast<seconds>(time_end - time_start).count() << "seconds" << endl;
}

vector<Ciphertext> expand_ciphertext(
    const Ciphertext &encrypted, uint32_t m, GaloisKeys &galkey, vector<int> &galois_elts)
{
    uint32_t logm = ceil(log2(m));
    Plaintext two("2");
    auto n = ckks->poly_modulus_degree();
    vector<Ciphertext> temp;
    temp.push_back(encrypted);
    Ciphertext tempctxt;
    Ciphertext tempctxt_rotated;
    Ciphertext tempctxt_shifted;
    Ciphertext tempctxt_rotatedshifted;

    for (uint32_t i = 0; i < logm - 1; i++) {
        vector<Ciphertext> newtemp(temp.size() << 1);
        int index_raw = (n << 1) - (1 << i);
        int index = (index_raw * galois_elts[i]) % (n << 1);
        for (uint32_t a = 0; a < temp.size(); a++) {
            ckks->evaluator->rotate_vector(temp[a], ckks->rots[i], ckks->galois_keys, tempctxt_rotated);
            ckks->evaluator->add(temp[a], tempctxt_rotated, newtemp[a]);
            multiply_power_of_X(temp[a], tempctxt_shifted, index_raw);
            multiply_power_of_X(tempctxt_rotated, tempctxt_rotatedshifted, index);
            ckks->evaluator->add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
        }
        temp = newtemp;
    }

    vector<Ciphertext> newtemp(temp.size() << 1);
    int index_raw = (n << 1) - (1 << (logm - 1));
    int index = (index_raw * galois_elts[logm - 1]) % (n << 1);
    for (uint32_t a = 0; a < temp.size(); a++) {
        if (a >= (m - (1 << (logm - 1)))) {
            ckks->evaluator->multiply_plain(temp[a], two, newtemp[a]);

        } else {
            ckks->evaluator->rotate_vector(temp[a], galois_elts[logm - 1], galkey, tempctxt_rotated);
            ckks->evaluator->add(temp[a], tempctxt_rotated, newtemp[a]);
            multiply_power_of_X(temp[a], tempctxt_shifted, index_raw);
            multiply_power_of_X(tempctxt_rotated, tempctxt_rotatedshifted, index);
            ckks->evaluator->add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
        }
    }

    vector<Ciphertext>::const_iterator first = newtemp.begin();
    vector<Ciphertext>::const_iterator last = newtemp.begin() + m;
    vector<Ciphertext> newVec(first, last);

    return newVec;
}

void multiply_power_of_X(const Ciphertext &encrypted, Ciphertext &destination, uint32_t index)
{
    auto coeff_mod_count = ckks->encryption_parameters->coeff_modulus().size() - 1;
    auto coeff_count = ckks->poly_modulus_degree();
    auto encrypted_count = encrypted.size();

    destination = encrypted;
    for (int i = 0; i < encrypted_count; i++) {
        for (int j = 0; j < coeff_mod_count; j++) {
            negacyclic_shift_poly_coeffmod(
                encrypted.data(i) + (j * coeff_count),
                coeff_count,
                index,
                ckks->get_modulus(destination, j),
                destination.data(i) + (j * coeff_count));
        }
    }
}

vector<double> gen_random_vector(size_t n, mt19937 &gen)
{
    uniform_int_distribution<> dis(1, 100);

    vector<double> vec(n);
    for (int i = 0; i < n; ++i) {
        vec[i] = dis(gen);
    }

    return vec;
}