#include "sgn.h"
#include <chrono>
#include <cassert>
#include <cstdio>
#include <mutex>
#include <iomanip>
#include <thread>

#include <seal/seal.h>


using namespace std;
using namespace seal;
// using namespace std::chrono;

// It's ok to have x as the dest, we just overwrite
void SgnEvaluator::sgn(int dg, int df, Ciphertext& x, Ciphertext& dest, shared_ptr<CKKSManager> ckks) {
    dest = x;
    for (int i = 1; i <= dg; i++) {
        if (i == dg && df == 0) {
            // Special coeffs for last poly 
            eval_odd_deg9_poly(g4_coeffs_last, dest, dest, ckks);
        } else {
            eval_odd_deg9_poly(g4_coeffs, dest, dest, ckks);
        }
    }
    for (int i = 1; i <= df; i++) {
        if (i == df) {
            // Special coeffs for last poly
            eval_odd_deg9_poly(f4_coeffs_last, dest, dest, ckks);
        } else {
            eval_odd_deg9_poly(f4_coeffs, dest, dest, ckks);
        }
    }
}



void SgnEvaluator::eval_odd_deg9_poly(vector<double>& a, Ciphertext& x, Ciphertext& dest, shared_ptr<CKKSManager> ckks) {
    bool DBG = false;
    /*
        (polyeval/odd9.h)
        P(x) = a9 x^9 + a7 x^7 + a5 x^5 + a3 x^3 + a1 x

        T1 = (a3 + a5 x^2) x^3 
        T2 = (a7 x + a9 x^3) x^6
        T3 = a1 x
        P(x) = T1 + T2 + T3

        Depth=4, #Muls=5

        Exactly what babystep_giantstep would do, but written explicitly to optimize 

        ###

        -> Errorless Polynomial Evaluation (3.2. of https://eprint.iacr.org/2020/1203)
        GOAL: evaluate a polynomial exactly so no need to stabilize and lose precision
        (x at level L and scale D --> P(x) at level L-4 and scale D)
        it's possible to do this exactly for polyeval as (x,x2,x3,x6) determine the scale D_L for each involved level L:
        (assume the primes at levels L to L-4 are p, q, r, s)

        level       ctx       scale (D_l)
        ==================================
          L          x          D
          L-1        x2         D^2 / p
          L-2        x3         D^3 / pq
          L-3        x6         D^6 / p^2 q^2 r

        Just by encoding constants at different scales we can make every ctx at level l be at scale D_l
        (not possible in general, e.g. rescale(x2*x2) produces L-2 ciphertext with scale D^4/ppq)
        (to fix this we would use the Adjust op. that multiplies ctx by constants and Algo 3 for primes from https://eprint.iacr.org/2020/1118)

        Now we know that sc(P(x)) should be D, so we recursively go back to compute the scales for each coefficient
        sc(T1)=sc(T2)=sc(T3)=sc(P(x))=D

        T3: 
            sc(a1) = q (should be p but it gets multiplied with modswitched x)

        T2:
            sc(x^6) = D^6 / p^2 q^2 r, so sc(a7*x) = sc(a9*x^3) = p^2 q^2 r s / D^5
            next, sc(a7) = p^2 q^3 r s / D^6
            similarly, sc(a9) = p^3 q^3 r^2 s / D^8
        
        T1:
            sc(x^3) = D^3 / pq
            implying sc(a3) = pqr / D^2 and also sc(a5*x^2) = pqr / D^2
            as sc(x^2) = D^2 / p this implies sc(a5) = p^2 q^2 r / D^4
    */
    // chrono::high_resolution_clock::time_point time_start, time_end;
    // time_start = high_resolution_clock::now();
    int n_ct_muls = 0;
    size_t level_start = x.coeff_modulus_size(); // L
    double D = x.scale(); // maybe not init_scale but preserved

    
    shared_ptr<Evaluator> evaluator = ckks->evaluator;

    uint64_t p = ckks->get_modulus(x, 1);
    uint64_t q = ckks->get_modulus(x, 2);
    uint64_t r = ckks->get_modulus(x, 3);
    uint64_t s = ckks->get_modulus(x, 4);
    uint64_t t = ckks->get_modulus(x, 5);

    p=q;
    q=r;
    r=s;
    s=t;

    vector<double> a_scales(10);
    a_scales[1] = q;
    a_scales[3] = (double)p/D * q/D * r;
    a_scales[5] = (double)p/D * p/D * q/D * q/D * r;
    a_scales[7] = (double)p/D * p/D * q/D * q/D * q/D * r/D * s;
    a_scales[9] = (double)p/D * p/D * p/D * q/D * q/D * q/D * r/D * r/D * s;

    ///////////////////////////////////////////////
    Ciphertext x2, x3, x6;

    evaluator->square(x, x2), n_ct_muls++;
    evaluator->relinearize_inplace(x2, ckks->relin_keys);
    evaluator->rescale_to_next_inplace(x2); // L-1

    evaluator->mod_switch_to_next_inplace(x); // L-1
    evaluator->multiply(x2, x, x3), n_ct_muls++;
    evaluator->relinearize_inplace(x3, ckks->relin_keys);  
    evaluator->rescale_to_next_inplace(x3); // L-2

    evaluator->square(x3, x6), n_ct_muls++;
    evaluator->relinearize_inplace(x6, ckks->relin_keys);
    evaluator->rescale_to_next_inplace(x6); // L-3

    if (DBG) cout << "[odd9] Constructed x^2, x^3, x^6" << endl;

    Plaintext a1, a3, a5, a7, a9;

    // Build T1
    Ciphertext T1;
    double a5_scale = D/x2.scale()*p/x3.scale()*q;
    ckks->encoder->encode(a[5], x2.parms_id(), a5_scale, a5); // L-1
    evaluator->multiply_plain(x2, a5, T1);
    evaluator->rescale_to_next_inplace(T1); // L-2

    // Update: using a_scales[3] is only approx. correct, so we directly use T1.scale()
    ckks->encoder->encode(a[3], T1.parms_id(), T1.scale(), a3); // L-2
    assert(fabs(T1.scale() - a3.scale()) < 1);
    
    evaluator->add_plain_inplace(T1, a3); // L-2

    evaluator->multiply_inplace(T1, x3), n_ct_muls++; 
    evaluator->relinearize_inplace(T1, ckks->relin_keys);
    evaluator->rescale_to_next_inplace(T1); // L-3


    // Build T2
    Ciphertext T2;
    Plaintext a9_switched;
    double a9_scale = D/x3.scale()*r/x6.scale()*q;
    ckks->encoder->encode(a[9], x3.parms_id(), a9_scale, a9); // L-2
    evaluator->multiply_plain(x3, a9, T2);
    evaluator->rescale_to_next_inplace(T2); // L-3

    Ciphertext a7x;
    double a7_scale = T2.scale()/x.scale()*p;
    ckks->encoder->encode(a[7], x.parms_id(), a7_scale, a7); // L-1 (x was modswitched)
    evaluator->multiply_plain(x, a7, a7x);
    evaluator->rescale_to_next_inplace(a7x); // L-2
    evaluator->mod_switch_to_inplace(a7x, T2.parms_id()); // L-3


    assert(fabs(T2.scale() - a7x.scale()) < 1);
    double mid_scale = (T2.scale() + a7x.scale()) / 2;
    T2.scale() = a7x.scale() = mid_scale; // this is the correct scale now, need to set it still to avoid SEAL assert
    evaluator->add_inplace(T2, a7x); // L-3
    evaluator->multiply_inplace(T2, x6), n_ct_muls++;
    evaluator->relinearize_inplace(T2, ckks->relin_keys);
    evaluator->rescale_to_next_inplace(T2); // L-4



    // Build T3
    Ciphertext T3;
    ckks->encoder->encode(a[1], x.parms_id(), p, a1); // L-1 (x was modswitched)
    evaluator->multiply_plain(x, a1, T3);
    evaluator->rescale_to_next_inplace(T3); // L-2

    // T1, T2 and T3 should be on the same scale up to floating point 
    // but we still need to set them manually to avoid SEAL assert
    assert(fabs(T1.scale() - D) < 1);
    assert(fabs(T2.scale() - D) < 1);
    assert(fabs(T3.scale() - D) < 1);
    double mid3_scale = (T1.scale() + T2.scale() + T3.scale()) / 3;
    T1.scale() = T2.scale() = T3.scale() = mid3_scale;

    dest = T2;
    evaluator->mod_switch_to_inplace(T1, dest.parms_id()); // L-4
    evaluator->add_inplace(dest, T1);
    evaluator->mod_switch_to_inplace(T3, dest.parms_id()); // L-4
    evaluator->add_inplace(dest, T3);

    /////////////////////////////////////////
    assert(level_start - dest.coeff_modulus_size() == 4);
    assert(n_ct_muls == 5);
    assert(fabs(dest.scale() - D) < 1);
    // it should be ==D but we don't stabilize if it's not, D' != D is ok
    // the goal was to make T1+T2+T3 work with minimal loss in precision
    // time_end = high_resolution_clock::now();
    // cout << "Poly eval took " << duration_cast<milliseconds>(time_end - time_start).count() << " ms" << endl;

}

void SgnEvaluator::eval_seg_gelu(vector<double>& a, Ciphertext& x, Ciphertext& dest_p, Ciphertext& dest_q, shared_ptr<CKKSManager> ckks) {

    shared_ptr<Evaluator> evaluator = ckks->evaluator;

    Ciphertext x2, x3, x4, x6;

    evaluator->square(x, x2);
    evaluator->relinearize_inplace(x2, ckks->relin_keys);
    evaluator->rescale_to_next_inplace(x2); // L-1

    evaluator->mod_switch_to_next_inplace(x); // L-1
    evaluator->multiply(x2, x, x3);
    evaluator->relinearize_inplace(x3, ckks->relin_keys);  
    evaluator->rescale_to_next_inplace(x3); // L-2

    evaluator->square(x2, x4);
    evaluator->relinearize_inplace(x4, ckks->relin_keys);
    evaluator->rescale_to_next_inplace(x4); // L-3

    evaluator->square(x3, x6);
    evaluator->relinearize_inplace(x6, ckks->relin_keys);
    evaluator->rescale_to_next_inplace(x6); // L-3

    // dest_p = a[3] * x^3 + a[2] * x^2 + a[1] * x + a[0]
    // dest_q = b[6] * x^6 + b[4] * x^4 + b[2] * x^2 + b[1] * x + b[0]
    vector<float> a_coeffs = {-0.5054031199708174, -0.42226581151983866, -0.11807612951181953, -0.011034134030615728};
    vector<float> b_coeffs = {0.008526321541038084,  0.5, 0.3603292692789629, 0.0, -0.037688200365904236, 0.0, 0.0018067462606141187};

    
}

