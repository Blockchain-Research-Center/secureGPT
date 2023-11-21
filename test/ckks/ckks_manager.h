#include <vector>
#include <seal/seal.h>


class CKKSManager {
public:
    CKKSManager(int N, int depth, double init_scale, bool VERBOSE) {
        this->N = N;
        this->depth = depth;
        this->init_scale = init_scale;
        init_ckks(VERBOSE);
    }

    uint64_t get_modulus(seal::Ciphertext& x, int k);
    void encode_and_encrypt(const std::vector<double>& msg, seal::Ciphertext& dest);
    void decrypt_and_decode(const seal::Ciphertext& ctx, std::vector<double>& dest);
    int poly_modulus_degree();
    int real_number_capacity();
    // Public for simplicity
    std::shared_ptr<seal::SEALContext> context;
    

    seal::PublicKey public_key;
    seal::SecretKey secret_key;
    seal::RelinKeys relin_keys;
    std::vector<int> rots;
    seal::GaloisKeys galois_keys;
    
    std::shared_ptr<seal::EncryptionParameters> encryption_parameters;
    std::shared_ptr<seal::CKKSEncoder> encoder;
    std::shared_ptr<seal::Encryptor> encryptor;
    std::shared_ptr<seal::Decryptor> decryptor;
    std::shared_ptr<seal::Evaluator> evaluator;

private:
    int N;
    int depth;
    double init_scale;

    void init_ckks(bool VERBOSE);
};