#pragma once

#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "globals.h"
#include <map>
#include <unordered_set>
#include <fstream> // Include the necessary header file


using namespace std;
using namespace sdsl;
// // Assumes that a root node always exists
// template <typename subset_rank_t>
//range struct
string bit_vector_to_string(const bit_vector& b){
    string s{};
    for(int i = 0; i < b.size(); i++){
        // cout<<b[i];
        if(b[i]) s.push_back('1');
        else s.push_back('0');
    }
    return s;
}
struct range{
    int64_t l;
    int64_t r;
    range(int64_t l, int64_t r): l(l), r(r){}
    range(){}
    bool operator==(const range& other) const{
        return l == other.l && r == other.r;
    }
    bool operator!=(const range& other) const{
        return !(*this == other);
    }
    int64_t size() const{
        return r-l+1;
    }
    bool empty() const{
        return r < l;
    }
    void set(int64_t l, int64_t r){
        this->l = l;
        this->r = r;
    }
    void set(const range& other){
        l = other.l;
        r = other.r;
    }
    void print() const{
        cout << "[" << l << ", " << r << "]";
    }
    string to_string() const{
        return "[" + std::to_string(l) + ", " + std::to_string(r) + "]";
    }
};

class SBWT_alphabet{

public:

//     subset_rank_t subset_rank; // The subset rank query implementation
    std::vector<int64_t> C; // The array of cumulative character counts
    std::vector<sd_vector<>> bit_vectors; // The bit vectors for each character
    std::vector<sd_vector<>::rank_1_type> rank_supports; // The rank support structures for each character
    std::vector<bit_vector> word_bit_vectors; // 1 if corresponding k-mer of that length is a word
    std::vector<sdsl::rank_support_v5<>> word_rank_supports;
    int k;
//     vector<pair<int64_t,int64_t> > kmer_prefix_precalc; // SBWT intervals for all p-mers with p = precalc_k.
//     int64_t precalc_k = 0;

//     int64_t n_nodes; // Number of nodes (= columns) in the data structure
//     int64_t n_kmers; // Number of k-mers indexed in the data structure
//     int64_t k; // The k-mer k

    int64_t get_char_idx(char c) const{
        if(!is_valid_char(c)) return -1;
        return c-'a';
    }
    bool is_valid_char(char c) const{
        return c >= 'a' && c <= 'z';
    }
    bool is_valid_string(const string& s) const{
        for(char c : s){
            if(!is_valid_char(c)) return false;
        }
        return true;
    }
    char get_char(int64_t idx) const{
        return 'a' + idx;
    }
    bool same_group(const string& s1, const string& s2, int k) const{
        if(s1.size()<k-1 || s2.size()<k-1) return false;
        for(int i = 0; i < k-1; i++){
            if(s1[s1.size()-i-1] != s2[s2.size()-i-1]) return false;
        }
        return true;
    }
    SBWT_alphabet(){
        C.resize(26+2);
        bit_vectors.resize(26);
        rank_supports.resize(26);
    }
    void build(const unordered_map<string,int64_t>& words, const vector<string>& sorted_kmers, const unordered_set<string>& kmers, int _k){
        k = _k;
        // Initialize the bit vectors
        auto __start = __get_now();
        int64_t n_kmers = sorted_kmers.size();
        for(int i = 0; i < 26; i++){
            if(i!=0) std::cerr << '\r' << std::flush; // Erase current line
            std::cerr << "bit vector " << get_char(i) << " ("<< i << "/26)" << std::flush;
            bit_vector b = bit_vector(n_kmers, 0);
            for(int j = 0; j < n_kmers; j++){
                if(!is_valid_string(sorted_kmers[j])) throw runtime_error("Invalid character in k-mer: " + sorted_kmers[j]);
                if((j==0)||(!same_group(sorted_kmers[j], sorted_kmers[j-1], k))){
                    //assume(d) no kmers are in the same k-1 suffix group
                    string new_str = sorted_kmers[j]+get_char(i);
                    if(new_str.size() > k) new_str = new_str.substr(1, k);
                    if(kmers.find(new_str) != kmers.end()){
                        b[j] = 1;
                    }
                }
                // else{
                //     cout << "same group: " << sorted_kmers[j] << " " << sorted_kmers[j-1] << endl;
                // }
            }
            bit_vectors[i] = sd_vector<>(b);
            rank_supports[i] = sd_vector<>::rank_1_type(&bit_vectors[i]);
        }
        std::cerr << '\r' << std::flush; // Erase current line
        auto __end = __get_now();
        __get_duration("bit vectors done", __start, __end);

        // Initialize the word bit vectors
        __start = __get_now();
        word_bit_vectors.resize(k);
        for(int i = 0; i < k; i++) word_bit_vectors[i] = bit_vector(n_kmers, 0);

        word_rank_supports.resize(k);
        for(int i = 0; i < n_kmers; i++){
            if(words.find(sorted_kmers[i]) != words.end()){
                word_bit_vectors[sorted_kmers[i].size()-1][i] = 1;
            }
        }
        for(int i = 0; i < k; i++){
            word_rank_supports[i] = rank_support_v5<>(&word_bit_vectors[i]);
        }
        __end = __get_now();
        __get_duration("word bit vectors done", __start, __end);

        // Compute the cumulative character counts
        C[0] = 0;
        C[1] = 1;//empty string
        for(int i = 0; i < 26; i++){
            C[i+2] = C[i+1] + rank_supports[i](n_kmers);
        }
        if(C[27] != n_kmers) throw runtime_error("Error: C[27] != n_kmers");
    }

    void save(const string& filename) const{
        std::ofstream out(filename, std::ios::binary); // Replace throwing_ofstream with std::ofstream
        for(int i = 0; i < 26; i++){
            bit_vectors[i].serialize(out);
        }
        out.write((char*)&k, sizeof(k));
        for(int i = 0; i < k; i++){
            word_bit_vectors[i].serialize(out);
        }
        //save k
    }
    void load(const string& filename){
        std::ifstream in(filename, std::ios::binary); // Replace throwing_ifstream with std::ifstream

        bit_vectors.resize(26);
        rank_supports.resize(26);

        for(int i = 0; i < 26; i++){
            bit_vectors[i].load(in);
            rank_supports[i] = sd_vector<>::rank_1_type(&bit_vectors[i]);
        }

        in.read((char*)&k, sizeof(k));
        cout << "k: " << k << endl;

        word_bit_vectors.resize(k);
        word_rank_supports.resize(k);
        for(int i = 0; i < k; i++){
            word_bit_vectors[i].load(in);
            word_rank_supports[i] = rank_support_v5<>(&word_bit_vectors[i]);
        }

        C.resize(26+2);
        C[0] = 0;
        C[1] = 1;
        for(int i = 0; i < 26; i++){
            C[i+2] = C[i+1] + rank_supports[i](bit_vectors[i].size());
        }

    }
    range forward(range I, char c) const{
        int64_t char_idx = get_char_idx(c);
        if(char_idx == -1) return range(-1,-2);
        // cout << "char_idx: " << char_idx << endl;
        // I= range(C[char_idx+1] + rank_supports[char_idx](I.l),
        // cout << "I: "; I.print(); cout << endl;
        // return I;
        return range(C[char_idx+1] + rank_supports[char_idx](I.l),
                    C[char_idx+1] + rank_supports[char_idx](I.r+1)-1);
    }
    int64_t search(const string& kmer) const{
        if(!is_valid_string(kmer)) return -1;
        range I = {0, C[27]-1};
        for(char c : kmer){
            I = forward(I, c);
            if(I.empty()) return -1;
        }
        return I.l;
    }
    // unordered_set<int64_t> anagram_search(const string& kmer, const unordered_map<string,int64_t>& words, ofstream& out_stream){
    unordered_set<int64_t> kangaroo_search(const string& kmer, const unordered_map<string,int64_t>& words)const{
        //clear positions, intervals and used_symbols
        vector<range> intervals{};
        vector<int> positions{};
        bit_vector used_symbols=bit_vector(kmer.size(), 0);

        unordered_set<int64_t> out{};
        
        if(!is_valid_string(kmer)) return {};
        
        string s{};
        range I = {0, C[27]-1};
        int i = 0;
        
        while(i < kmer.size() || s.size() > 0){
            if(i == kmer.size()){
                i = positions.back()+1;
                positions.pop_back();
                I = intervals.back();
                intervals.pop_back();
                s.pop_back();
                used_symbols[i-1] = 0;
                continue;
            }
            if(used_symbols[i]){
                i++;
                continue;
            }
            range I_new = forward(I, kmer[i]);
            if(!I_new.empty()){
                s.push_back(kmer[i]);
                if(word_rank_supports[s.size()-1](I_new.r+1) - word_rank_supports[s.size()-1](I_new.l) > 0){
                    int64_t pos = words.at(s);
                    if(out.find(pos) == out.end()){
                        out.insert(pos);
                    }else{
                        s.pop_back();
                        i++;
                        continue;
                    }
                }
                positions.push_back(i);
                intervals.push_back(I);
                I = I_new;
                used_symbols[i] = 1;
                i++;
            }else{
                i++;
            }
        }
        return out;
    }
    unordered_set<int64_t> anagram_search(const string& kmer, const unordered_map<string,int64_t>& words)const{
        if(!is_valid_string(kmer)) return {};
        vector<range> intervals{};
        vector<int> positions{};
        //get unique characters
        unordered_map<char,int> char_count;
        int64_t num_iterations = 0;
        int64_t average_depth = 0;
        int64_t num_attempts = 0;
        int64_t num_successes = 0;
        int64_t num_collisions = 0;
        for(char c : kmer){
            char_count[c]++;
        }
        vector<pair<char,int>> unique_chars;
        for(auto& p : char_count){
            unique_chars.push_back(p);
        }
        //sort by count
        sort(unique_chars.begin(), unique_chars.end(), [](const pair<char,int>& a, const pair<char,int>& b){
            // return a.second < b.second;
            return a.first < b.first;
        });
        unordered_set<int64_t> out{};
        string s{};
        range I = {0, C[27]-1};
        int i = 0;
        
        while(i < unique_chars.size() || s.size() > 0){
            num_iterations++;
            average_depth += s.size();
            if(i == unique_chars.size()){
                // printf("i: %d, s: %s, I: %s, bits: %s\n", i, s.c_str(), I.to_string().c_str(),bit_vector_to_string(used_symbols).c_str());
                i = positions.back()+1;
                positions.pop_back();
                I = intervals.back();
                intervals.pop_back();
                s.pop_back();
                unique_chars[i-1].second++;
                continue;
            }
            if(!unique_chars[i].second){
                i++;
                continue;
            }
            range I_new = forward(I, unique_chars[i].first);
            if(!I_new.empty()){
                // printf("i: %d, s: %s, I: %s, bits: %s\n", i, s.c_str(), I.to_string().c_str(), bit_vector_to_string(used_symbols).c_str());
                s.push_back(unique_chars[i].first);
                num_attempts++;
                if(word_rank_supports[s.size()-1](I_new.r+1) - word_rank_supports[s.size()-1](I_new.l) > 0){
                    // out.push_back(words.at(s));
                    int64_t pos = words.at(s);
                    num_successes++;
                    if(out.find(pos) == out.end()){
                        out.insert(pos);
                    }else{
                        num_collisions++;
                        s.pop_back();
                        i++;
                        continue;
                    }
                }
                positions.push_back(i);
                intervals.push_back(I);
                I = I_new;
                unique_chars[i].second--;
                i=0;
                // i++;
            }else{
                i++;
            }
        }
        // cout << "num_iterations: " << num_iterations << endl;
        // cout << "average_depth: " << (double)average_depth/num_iterations << endl;
        // cout << "num_attempts: " << num_attempts << endl;
        // cout << "num_successes: " << num_successes << endl;
        // cout << "num_collisions: " << num_collisions << endl;
        return out;
    }

};
