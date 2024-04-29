#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <sdsl/vectors.hpp>
#include <unordered_set>
#include <queue>
#include <omp.h>
#include <atomic>
#include <filesystem>
#include <thread>

// #include "Main/ColorSearchMain.h"
// #include "Main/IndexSearchMain.h"
// #include "Main/Main.h"
//!DEBUG
#include "globals.h"
#include "sbwt_alphabet.hpp"


// using sbwt_search::ColorSearchMain;
// using sbwt_search::IndexSearchMain;
// using sbwt_search::Main;
using std::cout;
// using gpu_utils::get_free_gpu_memory;
// using gpu_utils::get_taken_gpu_memory;
// using memory_utils::get_total_system_memory;
// using std::span; //!DEBUG
using std::string;
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include <random>
#define __size_mb(v) cout << #v << " size: " << size_in_mega_bytes(v) << "MB" << std::endl;



int main3()
{
    // auto end = __get_now();
    // __get_duration("Unique substrings", start, end);
    auto start = __get_now();
    bit_vector b = bit_vector(80*(1<<20), 0);
    for (size_t i=0; i < b.size(); i+=26)
        b[i] = 1;
    auto end = __get_now();
    __get_duration("make bit_vector", start, end);
    __size_mb(b);

    start = __get_now();
    sdsl::rank_support_v5<> b_rank(&b);
    end = __get_now();
    __get_duration("make rank_support_v5", start, end);
    __size_mb(b_rank);
    for (size_t i=0; i<=b.size(); i+= b.size()/4)
        cout << "(" << i << ", " << b_rank(i) << ") ";
    cout << endl;

    sd_vector<> sdb(b);
    __size_mb(sdb);

     start = __get_now();
    sd_vector<>::rank_1_type sdb_rank(&sdb);
     end = __get_now();
    __get_duration("make sd_vector", start, end);
    __size_mb(sdb_rank);
    for (size_t i=0; i<=sdb.size(); i+= sdb.size()/4)
        cout << "(" << i << ", " << sdb_rank(i) << ") ";
    cout << endl;

    rrr_vector<> rrrb(b);
    __size_mb(rrrb);

    start = __get_now();
    rrr_vector<>::rank_1_type rrrb_rank(&rrrb);
    end = __get_now();
    __get_duration("make rrr_vector", start, end);
    __size_mb(rrrb_rank);
    for (size_t i=0; i<=rrrb.size(); i+= rrrb.size()/4)
        cout << "(" << i << ", " << rrrb_rank(i) << ") ";
    cout << endl;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, b.size() - 1);
    const int num_queries = 10000000;
     start = __get_now();
    for (int i=0; i<num_queries; i++)
        b_rank(dis(gen));
     end = __get_now();
    __get_duration("Rank queries on bit_vector", start, end);

    start = __get_now();
    for (int i=0; i<num_queries; i++)
        sdb_rank(dis(gen));
    end = __get_now();
    __get_duration("Rank queries on sd_vector", start, end);

    start = __get_now();
    for (int i=0; i<num_queries; i++)
        rrrb_rank(dis(gen));
    end = __get_now();
    __get_duration("Rank queries on rrr_vector", start, end);
    return 0;
}

enum class leaderboard_type{
    anagram,
    kangaroo,
    n_kangaroo
};

bool global_filter_condition(const std::string& str) {
    // return str.size() > 3;
    return true;
}
bool filter_condition(const std::string& str) {
    return str.size() > 2;
    // return true;
}
double weight_condition(const std::string& str, const std::string& parent) {
    return ((double)str.size())*str.size()/parent.size();
    // return ((double)str.size())/parent.size();
// if(subword in word) return 0
// else return len(subword)/len(word)
    // if(parent.find(str)!=std::string::npos) return 0;//when subword in word
    // return ((double)str.size())/parent.size();
}

struct PriorityQueueItem {
    double weighted_size;
    size_t size;
    std::string str;
    bool operator<(const PriorityQueueItem& other) const {
        // return size > other.size;
        return weighted_size > other.weighted_size;
    }
};
std::mutex cout_mutex;
std::atomic<int> progress(0);
std::atomic<int64_t> total_subword(0);
std::atomic<int> kangaroo_num(2);
void processLines(const std::vector<std::string>& lines, 
                  const SBWT_alphabet& sbwt, 
                  const std::unordered_map<string,int64_t>& unique_words,
                  std::vector<std::priority_queue<PriorityQueueItem>>& pqs,
                  size_t start, size_t end, int thread_id, 
                  int64_t global_time,
                //   bool do_kangaroo_search = 1) {
                  leaderboard_type lb_type = leaderboard_type::kangaroo) {
    for(size_t i = start; i < end; i++) {
        // auto matches = sbwt.kangaroo_search(lines[i], unique_words);
        unordered_set<int64_t> matches;
        unordered_set<MatchCollection,MatchCollectionHash> matches2;
        if(lb_type==leaderboard_type::kangaroo){
            matches = sbwt.kangaroo_search(lines[i], unique_words);
        } else if(lb_type==leaderboard_type::n_kangaroo){
            matches2 = sbwt.n_kangaroo_search(lines[i], unique_words, kangaroo_num);
        }else{
            matches = sbwt.anagram_search(lines[i], unique_words);
        }
        // auto matches = sbwt.anagram_search(lines[i], unique_words);
        // size_t true_count = matches.size();
        size_t true_count = 0;
        double weighted_count = 0;
        if(!(lb_type==leaderboard_type::n_kangaroo)){
            total_subword += matches.size();
             for(auto match : matches) true_count+=filter_condition(lines[match]);
             for(auto match : matches) weighted_count+=filter_condition(lines[match])*weight_condition(lines[match], lines[i]);
        }else{
            true_count+=matches2.size();
            for(auto coll : matches2){
                for(int64_t m : coll.v){
                    weighted_count+=weight_condition(lines[m], lines[i]);
                }
            }
            // if(lines[i]=="cholestyramines"){
            //     cout << "cholestyramines: " << true_count << " " << weighted_count << std::endl;
            //     for(auto coll : matches2){
            //         for(int64_t m : coll.v){
            //             cout << lines[m] << " ";
            //         }
            //         cout << std::endl;
            //     }
            // }
        }
        if(pqs[thread_id].size() < 200) {
            pqs[thread_id].push({weighted_count,true_count, lines[i]});
        } else if(weighted_count > pqs[thread_id].top().weighted_size) {
            pqs[thread_id].pop();
            pqs[thread_id].push({weighted_count,true_count, lines[i]});
        }
        int p = ++progress;
        if (p % 500 == 0 && p != 0) {
            auto now = __get_now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count() - global_time;
            int64_t expected_duration = duration * lines.size() / p;
            int64_t expected_seconds_left = (expected_duration - duration) / 1000;
            std::lock_guard<std::mutex> guard(cout_mutex);
            std::cerr << "\rThread " << thread_id << ": Line " << p << "/" << lines.size() << ", ETA: " << expected_seconds_left << "s" << std::flush;
        }
    }
}


priority_queue<PriorityQueueItem> calc_leaderboard(const std::vector<std::string>& lines, 
                  const SBWT_alphabet& sbwt, 
                  const std::unordered_map<string,int64_t>& unique_words,
                  int num_threads,
                //   bool do_kangaroo_search = 1) {
                  leaderboard_type lb_type = leaderboard_type::kangaroo) {
        progress = 0;
        total_subword = 0;
            std::vector<std::thread> threads(num_threads);
        std::vector<std::priority_queue<PriorityQueueItem>> pqs(num_threads);
        auto start = __get_now();
        int64_t global_time = std::chrono::duration_cast<std::chrono::milliseconds>(start.time_since_epoch()).count();
        size_t lines_per_thread = lines.size() / num_threads;
        for(size_t t = 0; t < num_threads; t++) {
            size_t start = t * lines_per_thread;
            size_t end = (t == num_threads - 1) ? lines.size() : (t + 1) * lines_per_thread;
            threads[t] = std::thread(processLines, std::ref(lines), std::ref(sbwt), std::ref(unique_words), std::ref(pqs), start, end, t, global_time, lb_type);
        }
        for(auto& thread : threads) {
            thread.join();
        }
        cout << std::endl;
        priority_queue<PriorityQueueItem> pq;
        for(auto& pq_thread : pqs){
            while(!pq_thread.empty()){
                if(pq.size()<200){
                    pq.push(pq_thread.top());
                }else if(pq_thread.top().weighted_size>pq.top().weighted_size){
                    pq.pop();
                    pq.push(pq_thread.top());
                }
                pq_thread.pop();
            }
        }
    return pq;
}

bool is_in_alphabet(char c) {
    return (c >= 'a' && c <= 'z');
}

int main2(int argc, char **argv) {
    if(argc < 2) {
        std::cout << "Please provide a file name." << std::endl;
        return 1;
    }

    std::ifstream file(argv[1]);
    std::filesystem::path p(argv[1]);
    if(!file) {
        std::cout << "Could not open file." << std::endl;
        return 1;
    }

    std::vector<std::string> lines;
    std::string line;
    while(std::getline(file, line)) {
        // lines.push_back(line);
        //lowercase the line
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);
    if(global_filter_condition(line)){
            lines.push_back(line);
        }
    }

    std::set<char> unique_chars;
    for(const auto& str : lines) {
        unique_chars.insert(str.begin(), str.end());
    }


    int n=0;
    int m=0;
    for(char c : unique_chars){
        if(is_in_alphabet(c)) n++;
        else m++;
    }
    if(n==26&&m==0){
        std::cout << "Alphabet is complete" << std::endl;
    }else{
        std::cout << "Alphabet is invalid:" << std::endl;
        for(char c : unique_chars) {
            std::cout << c << ' ';
        }
        std::cout << std::endl;
        std::cout << "Skipping invalid strings" << std::endl;
        unique_chars.clear();
        int64_t skip_count = 0;
        std::vector<std::string> new_lines;
        for(const auto& str : lines) {
            bool valid = true;
            for(char c : str) {
                if(!is_in_alphabet(c)) {
                    valid = false;
                    skip_count++;
                    break;
                }
            }
            if(valid) {
                unique_chars.insert(str.begin(), str.end());
                new_lines.push_back(str);
            }
        }
        std::cout << "Skipped " << skip_count << " strings" << std::endl;
        lines = new_lines;
        std::cout << "New number of strings: " << lines.size() << std::endl;
        if(unique_chars.size() == 26) {
            std::cout << "Alphabet is now complete" << std::endl;
        } else {
            std::cout << "Alphabet is still invalid, continuing anyway" << std::endl;
        }
    }

    auto longest = std::max_element(lines.begin(), lines.end(), [](const std::string& a, const std::string& b) {
        return a.size() < b.size();
    });
    std::cout << "Longest string: " << *longest << std::endl;
    int k=longest->size();

    std::cout << std::endl;
    std::vector<std::vector<std::string>> bins(k + 1);
    for(int64_t i = 0; i < lines.size(); i++) {
        bins[lines[i].size()].push_back(lines[i]);
    }
    int64_t total = 0;
    for(size_t i = 0; i < bins.size(); ++i) {
        // std::cout << "Bin " << i << ": " << bins[i].size() << " strings" << std::endl;
        total += bins[i].size();
    }
    std::cout << "Total number of strings: " << total << std::endl;
    auto start = __get_now();
    std::unordered_set<std::string> unique_substrings;
    std::unordered_map<std::string,int64_t> unique_words;
    unique_substrings.insert("");
    // for(const auto& str : lines) {
    for(int64_t i = 0; i < lines.size(); i++) {
        for(size_t j = 1; j <= lines[i].size(); ++j) {
            unique_substrings.insert(lines[i].substr(0, j));
        }
        // unique_words.insert(str);
        unique_words[lines[i]]=i;
    }
    auto end = __get_now();
    __get_duration("Unique words and substrings", start, end);

    std::cout << "Number of unique substrings: " << unique_substrings.size() << std::endl;
    start = __get_now();
    // sort colexicographically
    std::vector<std::string> sorted_substrings(unique_substrings.begin(), unique_substrings.end());

    for(auto& str : sorted_substrings) std::reverse(str.begin(), str.end());
    std::sort(sorted_substrings.begin(), sorted_substrings.end());
    for(auto& str : sorted_substrings) std::reverse(str.begin(), str.end());
    end = __get_now();

    __get_duration("Sort substrings", start, end);
    if(sorted_substrings.size()!=unique_substrings.size()){
        std::cout << "Error: sorted_substrings.size()!=unique_substrings.size()" << std::endl;
    }
    SBWT_alphabet sbwt=SBWT_alphabet();


    // string sbwt_file = "C:\\Code\\SBWT_alphabet\\sbwt_alphabet2.sdsl";
    std::string sbwt_file = p.replace_extension(".sbwt_alphabet.sdsl").string();
    //if does not exist, build it, otherwise load it
    if(!std::filesystem::exists(sbwt_file)){
        cout << "Building SBWT to " << sbwt_file << std::endl;
        start = __get_now();
        sbwt.build(unique_words, sorted_substrings, unique_substrings, k);
        sbwt.save(sbwt_file);
        end = __get_now();
        __get_duration("Build SBWT", start, end);
    }else{
        cout << "Loading SBWT from " << sbwt_file << std::endl;
        start = __get_now();
        sbwt.load(sbwt_file);
        end = __get_now();
        __get_duration("Load SBWT", start, end);
    }
    // sbwt.build(unique_words, sorted_substrings, unique_substrings, k);
    // sbwt.save(sbwt_file);

    // start = __get_now();
    // sbwt.load(sbwt_file);
    // end = __get_now();
    // __get_duration("Load SBWT", start, end);

    // //sample some strings
    // for(int i=0; i<10; i++){
    // // for(int i=0; i<1; i++){
    //     string s=lines[i];
    //     // string s="aaa";
    //     int64_t match=sbwt.search(s);
    //     if(match==-1){
    //         std::cout<<s<<": not found"<<std::endl;
    //     }else{
    //         cout<<s<<": "<<match<<" = "<<sorted_substrings[match]<<std::endl;
    //     }
    // }
    //test all strings
    start = __get_now();
    for(int i=0; i<lines.size(); i++){
    // for(int i=2000; i<2650; i++){
        if (i % 10000 == 0) {
            if(i!=0) std::cerr << '\r' << std::flush;
            std::cerr << "Testing " << i << "/" << lines.size() << std::flush;
        }
        string s=lines[i];
        // std::cerr << i << " " << s   << std::endl << std::flush;
        int64_t match=sbwt.search(s);
        if(match==-1){
            std::cout<<s<<": not found"<<std::endl;
            exit(1);
            break;
        }
        if(!sbwt.word_bit_vectors[s.size()-1][match]){
            std::cout<<s<<": bitvector not set"<<std::endl;
            exit(1);
            break;
        }
    }
    std::cerr << std::endl;
    end = __get_now();
    __get_duration("Search all strings", start, end);
    //print out the first couple of strings
    // for(int i=0; i<10; i++){
    //     cout<<sorted_substrings[i]<<std::endl;
    // }
    //whether aaa is in sorted_substrings
    // auto it=std::lower_bound(sorted_substrings.begin(), sorted_substrings.end(), "aaa");
    // if(it!=sorted_substrings.end()&&*it=="aaa"){
    //     cout<<"aaa is in sorted_substrings"<<std::endl;
    // }else{
    //     cout<<"aaa is not in sorted_substrings"<<std::endl;
    // }
    // vector<int64_t> matches = sbwt.anagram_search("chicken", unique_words);
    // for(int64_t match : matches){
    //     cout<<match<<std::endl;
    // }
    // for(int64_t match : matches){
    //     cout<<lines[match]<<std::endl;
    // }
    auto comp = [](pair<int64_t, string> a, pair<int64_t, string> b){
        return a.first>b.first;
    };
    //=
    priority_queue<pair<int64_t, string>, vector<pair<int64_t, string>>, decltype(comp)> pq(comp);
    start = __get_now();
    // for(int i=0; i<lines.size(); i++){
    // for(int i=0; i<1000; i++){
        // string &str=lines[i];
        // if (i % 500 == 0&&i!=0) {
        //     if(i!=0) std::cerr << '\r' << std::flush;
        //     auto now = __get_now();
        //     auto duration = chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
        //     auto expected_duration = duration * lines.size() / i;
        //     auto expected_seconds_left = (expected_duration - duration) / 1000;
        //     std::cerr << "String " << i << "/" << lines.size() << ", ETA: " << expected_seconds_left << "s" << std::flush;
        // }
    bool do_single_word_search = 0;
    bool do_leaderboard = 1;
    bool do_kangaroo_leaderboard = 0;
    bool do_anagram_leaderboard = 0;
    kangaroo_num=2;
    bool do_n_kangaroo_leaderboard = 1;
    if (do_single_word_search) {
        string str="softheartedness";
        if(argc==3){
            str=argv[2];
        }

        bool do_kangaroo_search = 0;
        bool do_anagram_search = 0;
        bool do_n_kangaroo_search = 1;
        // unordered_set<int64_t> matches = sbwt.anagram_search(str, unique_words);

        //=
        // store to file
        if(do_kangaroo_search){
            unordered_set<int64_t> matches = sbwt.kangaroo_search(str, unique_words);
            // ofstream out_stream("C:\\Code\\SBWT_alphabet\\anagrams4.txt");
            ofstream out_stream((p.parent_path() / "kangaroo.txt").string());
            for(int64_t match : matches){
                if(filter_condition(lines[match]))
                    out_stream<<lines[match]<<std::endl;
            }
            out_stream.close();
        }
        if(do_anagram_search){
            unordered_set<int64_t> matches = sbwt.anagram_search(str, unique_words);
            // ofstream out_stream("C:\\Code\\SBWT_alphabet\\anagrams5.txt");
            ofstream out_stream((p.parent_path() / "anagrams.txt").string());
            for(int64_t match : matches){
                if(filter_condition(lines[match]))
                    out_stream<<lines[match]<<std::endl;
            }
            out_stream.close();
        }

        if(do_n_kangaroo_search){
            unordered_set<MatchCollection,MatchCollectionHash> matches = sbwt.n_kangaroo_search(str, unique_words, kangaroo_num);
            // ofstream out_stream("C:\\Code\\SBWT_alphabet\\anagrams4.txt");
            ofstream out_stream((p.parent_path() / "kangaroo_covers.txt").string());
            for(auto coll : matches){
                for(int64_t m : coll.v){
                    out_stream<<lines[m] << " ";
                }
                out_stream<<std::endl;
            }
            out_stream.close();
        }
        //=

        // if(matches.size()>max_len_match){
        //     max_len_match=matches.size();
        //     max_len_str=str;
        // }
        // if(matches.size()>pq.top().first){
        // if(pq.size()<200){
        //     pq.push({matches.size(), str});
        // }else if(matches.size()>pq.top().first){
        //     pq.pop();
        //     pq.push({matches.size(), str});
        // }
        end = __get_now();
        __get_duration("Single word search", start, end);
    }
    //=

    if(do_leaderboard){
        int num_threads = 8;
        // int num_threads = std::thread::hardware_concurrency();
        cout << "Number of threads: " << num_threads << std::endl;
        if(do_kangaroo_leaderboard){
            priority_queue<PriorityQueueItem> pq=calc_leaderboard(lines, sbwt, unique_words, num_threads, leaderboard_type::kangaroo);
            // ofstream out_stream("C:\\Code\\SBWT_alphabet\\anagrams6.txt");
            ofstream out_stream((p.parent_path() / "kangaroo_leaderboard.txt").string());
            while(!pq.empty()){
                // out_stream<<pq.top().str<<": "<<pq.top().size<<std::endl;
                out_stream<<pq.top().str<<": "<<pq.top().size<<", ("<<pq.top().weighted_size<<")"<<std::endl;
                pq.pop();
            }
            cout<<"Total kangaroo words: "<<total_subword<<std::endl;
        }
        if(do_anagram_leaderboard){
            priority_queue<PriorityQueueItem> pq=calc_leaderboard(lines, sbwt, unique_words, num_threads, leaderboard_type::anagram);
            // ofstream out_stream("C:\\Code\\SBWT_alphabet\\anagrams7.txt");
            ofstream out_stream((p.parent_path() / "anagrams_leaderboard.txt").string());
            while(!pq.empty()){
                // out_stream<<pq.top().str<<": "<<pq.top().size<<std::endl;
                out_stream<<pq.top().str<<": "<<pq.top().size<<", ("<<pq.top().weighted_size<<")"<<std::endl;
                pq.pop();
            }
            cout<<"Total sub-anagrams: "<<total_subword<<std::endl;
        }
        if(do_n_kangaroo_leaderboard){
            priority_queue<PriorityQueueItem> pq=calc_leaderboard(lines, sbwt, unique_words, num_threads, leaderboard_type::n_kangaroo);
            // ofstream out_stream("C:\\Code\\SBWT_alphabet\\anagrams7.txt");
            ofstream out_stream((p.parent_path() / "n_kangaroo_leaderboard.txt").string());
            while(!pq.empty()){
                // out_stream<<pq.top().str<<": "<<pq.top().size<<std::endl;
                out_stream<<pq.top().str<<": "<<pq.top().size<<", ("<<pq.top().weighted_size<<")"<<std::endl;
                pq.pop();
            }
            cout<<"Total kangaroo covers: "<<total_subword<<std::endl;
        }

        end = __get_now();
        __get_duration("Kangaroo/Anagram search", start, end);
        // cout<<max_len_str<<": "<<max_len_match<<std::endl;
    }
    return 0;
}


int main(int argc, char **argv) {
    if(1){
        return main2(argc, argv);
    }else{
        return main3();
    }
    return 0;
}


//cmake ..
//cmake --build .
//.\bin\Debug\sbwt.exe ..\words_alpha.txt
//cmake .. ; cmake --build . ; .\bin\Debug\sbwt.exe ..\words_alpha.txt