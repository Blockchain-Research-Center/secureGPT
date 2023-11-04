#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <seal/seal.h>
#include <seal/util/polyarithsmallmod.h>

using namespace std::chrono;
using namespace seal::util;
using namespace std;
using namespace seal;

size_t poly_modulus_degree;
EncryptionParameters* params;
Evaluator* evaluator;
chrono::high_resolution_clock::time_point time_start, time_end;

vector<uint64_t> gen_random_vector(size_t n, mt19937 &gen);
vector<Ciphertext> expand_ciphertext(const Ciphertext &encrypted, uint32_t m, GaloisKeys &galkey, vector<uint32_t> galois_elts);
void multiply_power_of_X(const Ciphertext &encrypted, Ciphertext &destination, uint32_t index);

int main() {
    mt19937 gen(random_device{}());
    size_t poly_modulus_degree = 4096;
    params = new EncryptionParameters(scheme_type::bgv);
    params->set_poly_modulus_degree(poly_modulus_degree);
    params->set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    params->set_plain_modulus(786433);
    KeyGenerator keygen(*params);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    GaloisKeys gal_keys;
    vector<uint32_t> galois_elts;
    for (int i = 0; i < ceil(log2(poly_modulus_degree)); i++) {
        galois_elts.push_back((poly_modulus_degree + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
    }
    keygen.create_galois_keys(galois_elts, gal_keys);
    Encryptor encryptor(*params, public_key);
    Decryptor decryptor(*params, secret_key);
    evaluator = new Evaluator(*params);

    vector<Plaintext> a_pts;
    a_pts.reserve(768);
    for(int i = 0; i < 768; i++) {
        a_pts.emplace_back(gen_random_vector(poly_modulus_degree, gen));
    }

    vector<Ciphertext> b_compressed_cts;
    for(int i = 0; i < 768*768/poly_modulus_degree; i++) {
        Plaintext pt(gen_random_vector(poly_modulus_degree, gen));
        Ciphertext ct;
        encryptor.encrypt(pt, ct);
        b_compressed_cts.push_back(ct);
    }

    time_start = high_resolution_clock::now();
    vector<Ciphertext> b_expanded_cts;
    for (Ciphertext &ct : b_compressed_cts) {
        vector<Ciphertext> temp_cts = expand_ciphertext(ct, poly_modulus_degree, gal_keys, galois_elts);
        b_expanded_cts.insert(b_expanded_cts.end(), 
                              make_move_iterator(temp_cts.begin()), 
                              make_move_iterator(temp_cts.end()));
    }
    time_end = high_resolution_clock::now();
    cout << "expanding time: " << duration_cast<seconds>(time_end - time_start).count() << "seconds" << endl;

    vector<Ciphertext> res_cts;
    Plaintext pt(poly_modulus_degree, 0);
    time_start = high_resolution_clock::now();
    for(int i = 0; i < 768; i++) {
      Ciphertext res_col_ct;
      encryptor.encrypt(pt, res_col_ct);
      for(int j = 0; j < 768; j++) {
        Ciphertext temp;
        evaluator->multiply_plain(b_expanded_cts[i*768+j], a_pts[j], temp);
        evaluator->add(res_col_ct, temp, res_col_ct);
      }
      res_cts.push_back(res_col_ct);
    }
    time_end = high_resolution_clock::now();
    cout << "calculating res time: " << duration_cast<seconds>(time_end - time_start).count() << "seconds" << endl;
}

vector<uint64_t> gen_random_vector(size_t n, mt19937 &gen) {
    uniform_int_distribution<> dis(1, 100);  

    vector<uint64_t> vec(n);
    for (int i = 0; i < n; ++i) {
        vec[i] = dis(gen); 
    }
    
    return vec;
}

vector<Ciphertext> expand_ciphertext(const Ciphertext &encrypted, uint32_t m, GaloisKeys &galkey, vector<uint32_t> galois_elts) {

  uint32_t logm = ceil(log2(m));
  Plaintext two("2");
  auto n = params->poly_modulus_degree();
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
      evaluator->apply_galois(temp[a], galois_elts[i], galkey,tempctxt_rotated);
      evaluator->add(temp[a], tempctxt_rotated, newtemp[a]);
      multiply_power_of_X(temp[a], tempctxt_shifted, index_raw);
      multiply_power_of_X(tempctxt_rotated, tempctxt_rotatedshifted, index);
      evaluator->add(tempctxt_shifted, tempctxt_rotatedshifted,
                      newtemp[a + temp.size()]);
    }
    temp = newtemp;
  }

  vector<Ciphertext> newtemp(temp.size() << 1);
  int index_raw = (n << 1) - (1 << (logm - 1));
  int index = (index_raw * galois_elts[logm - 1]) % (n << 1);
  for (uint32_t a = 0; a < temp.size(); a++) {
    if (a >= (m - (1 << (logm - 1)))) {
      evaluator->multiply_plain(temp[a], two,
                                 newtemp[a]);

    } else {
      evaluator->apply_galois(temp[a], galois_elts[logm - 1], galkey,
                               tempctxt_rotated);
      evaluator->add(temp[a], tempctxt_rotated, newtemp[a]);
      multiply_power_of_X(temp[a], tempctxt_shifted, index_raw);
      multiply_power_of_X(tempctxt_rotated, tempctxt_rotatedshifted, index);
      evaluator->add(tempctxt_shifted, tempctxt_rotatedshifted,
                      newtemp[a + temp.size()]);
    }
  }

  vector<Ciphertext>::const_iterator first = newtemp.begin();
  vector<Ciphertext>::const_iterator last = newtemp.begin() + m;
  vector<Ciphertext> newVec(first, last);

  return newVec;
}

void multiply_power_of_X(const Ciphertext &encrypted, Ciphertext &destination, uint32_t index) {

  auto coeff_mod_count = params->coeff_modulus().size() - 1;
  auto coeff_count = params->poly_modulus_degree();
  auto encrypted_count = encrypted.size();

  destination = encrypted;
  for (int i = 0; i < encrypted_count; i++) {
    for (int j = 0; j < coeff_mod_count; j++) {
      negacyclic_shift_poly_coeffmod(encrypted.data(i) + (j * coeff_count),
                                     coeff_count, index,
                                     params->coeff_modulus()[j],
                                     destination.data(i) + (j * coeff_count));
    }
  }

}