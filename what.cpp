    int num_threads = 6;
    cout << "Number of threads: " << num_threads << std::endl;
    std::atomic<int> progress(0);
    vector<priority_queue<pair<int64_t, string>, vector<pair<int64_t, string>>, decltype(comp)>> pqs(num_threads, 
    priority_queue<pair<int64_t, string>, vector<pair<int64_t, string>>, decltype(comp)>(comp));

    #pragma omp parallel for
    for(int i=0; i<lines.size(); i++){
        string str=lines[i];
        unordered_set<int64_t> matches = sbwt.kangaroo_search(str, unique_words);

        int thread_id = omp_get_thread_num();
        if(pqs[thread_id].size()<200){
            pqs[thread_id].push({matches.size(), str});
        }else if(matches.size()>pqs[thread_id].top().first){
            pqs[thread_id].pop();
            pqs[thread_id].push({matches.size(), str});
        }
        if (i % 500 == 0 && i != 0) {
            auto now = __get_now();
            auto duration = chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
            auto expected_duration = duration * lines.size() / i;
            auto expected_seconds_left = (expected_duration - duration) / 1000;
            if(i!=0) std::cerr << '\r' << std::flush;
            std::cerr << "String " << i << "/" << lines.size() << ", ETA: " << expected_seconds_left << "s" << std::flush;
        }
    }
    
    cout << "Number of threads: " << num_threads << std::endl;
    priority_queue<pair<int64_t, string>, vector<pair<int64_t, string>>, decltype(comp)> pq(comp);
    for(auto& pq_thread : pqs){
        while(!pq_thread.empty()){
            if(pq.size()<200){
                pq.push(pq_thread.top());
            }else if(pq_thread.top().first>pq.top().first){
                pq.pop();
                pq.push(pq_thread.top());
            }
            pq_thread.pop();
        }
    }






    unordered_set<int64_t> kangaroo_search(const string& kmer, const unordered_map<string,int64_t>& words)const{
        //clear positions, intervals and used_symbols
        vector<range> intervals{};
        vector<int> positions{};
        bit_vector used_symbols=bit_vector(kmer.size(), 0);

        int64_t num_iterations = 0;
        int64_t average_depth = 0;
        int64_t num_attempts = 0;
        int64_t num_successes = 0;
        int64_t num_collisions = 0;
        // vector<int64_t> out{};
        unordered_set<int64_t> out{};
        
        if(!is_valid_string(kmer)) return {};
        
        string s{};
        range I = {0, C[27]-1};
        int i = 0;
        
        while(i < kmer.size() || s.size() > 0){
            num_iterations++;
            average_depth += s.size();
            if(i == kmer.size()){
                // printf("i: %d, s: %s, I: %s, bits: %s\n", i, s.c_str(), I.to_string().c_str(),bit_vector_to_string(used_symbols).c_str());
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
                num_attempts++;
                // printf("i: %d, s: %s, I: %s, bits: %s\n", i, s.c_str(), I.to_string().c_str(), bit_vector_to_string(used_symbols).c_str());
                s.push_back(kmer[i]);
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
                used_symbols[i] = 1;
                // i=0;
                i++;
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