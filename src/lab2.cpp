#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <cmath>

using namespace std;

// 动态k-mer选择函数，改进为更平滑的增长
int selectKmerSize(size_t ref_length) {
    if (ref_length < 500) return 10;
    if (ref_length < 2000) return 12;
    if (ref_length < 10000) return 14;
    if (ref_length < 50000) return 16;
    if (ref_length < 200000) return 18;
    return 20;
}

// 全局参数
int K;
const int MAX_EXTEND = 200; // 增加扩展长度以捕获更多匹配
const int MIN_ANCHOR_LENGTH = 15; // 降低最小锚点长度以包含更多潜在匹配
const int MISMATCH_TOLERANCE = 3; // 增加错配容忍度
const double SIMILARITY_THRESHOLD = 0.5; // 降低相似度阈值以提高灵敏度
const int GAP_PENALTY_BASE = 5; // 基础间隔惩罚
const double GAP_PENALTY_SCALE = 0.1; // 间隔惩罚缩放因子

// 获取反向互补序列
string getRC(const string &s) {
    string rc;
    rc.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        switch (*it) {
            case 'A': rc += 'T'; break;
            case 'T': rc += 'A'; break;
            case 'C': rc += 'G'; break;
            case 'G': rc += 'C'; break;
            default: rc += 'N'; break; // 处理非标准碱基
        }
    }
    return rc;
}

struct Anchor {
    int q_start;
    int q_end;
    int r_start;
    int r_end;
    int score;
    bool is_rc;
    
    Anchor(int qs, int qe, int rs, int re, bool rc = false) 
        : q_start(qs), q_end(qe), r_start(rs), r_end(re), is_rc(rc),
          score(calculateScore(qs, qe, rs, re)) {}
    
    bool operator<(const Anchor &other) const {
        return q_start < other.q_start;
    }
    
    bool overlaps(const Anchor &other) const {
        int q_overlap = min(q_end, other.q_end) - max(q_start, other.q_start);
        int r_overlap = min(r_end, other.r_end) - max(r_start, other.r_start);
        return q_overlap > 0 && r_overlap > 0;
    }
    
    // 计算相似度
    double similarity(const string &query, const string &ref) const {
        int matches = 0;
        int len = q_end - q_start;
        for (int i = 0; i < len; ++i) {
            char q_base = query[q_start + i];
            char r_base = is_rc ? getRC(string(1, ref[r_start + i]))[0] : ref[r_start + i];
            if (q_base == r_base) matches++;
        }
        return static_cast<double>(matches) / len;
    }
    
public:
    static int calculateScore(int qs, int qe, int rs, int re) {
        int len = min(qe - qs, re - rs);
        // 改进分数计算：对长匹配给予更高权重，同时考虑对齐长度
        return static_cast<int>(len * (1 + log2(1 + len)));
    }
};

// 并行构建k-mer索引
void buildKmerIndex(const string &seq, unordered_map<string, vector<int>> &kmerMap, int k, int start, int end, mutex &mtx) {
    unordered_map<string, vector<int>> localMap;
    for (int i = start; i + k <= end && i + k <= seq.size(); ++i) {
        string kmer = seq.substr(i, k);
        localMap[kmer].push_back(i);
    }
    lock_guard<mutex> lock(mtx);
    for (auto &pair : localMap) {
        kmerMap[pair.first].insert(kmerMap[pair.first].end(), pair.second.begin(), pair.second.end());
    }
}

vector<Anchor> findAnchors(const string &query, const string &ref) {
    vector<Anchor> anchors;
    
    // 动态选择k-mer大小
    K = selectKmerSize(ref.size());
    
    // 并行构建k-mer索引
    unordered_map<string, vector<int>> kmerMap;
    kmerMap.reserve(ref.size() - K + 1);
    mutex mtx;
    
    const int num_threads = thread::hardware_concurrency();
    vector<thread> threads;
    int chunk_size = ref.size() / num_threads;
    
    for (int i = 0; i < num_threads; ++i) {
        int start = i * chunk_size;
        int end = (i == num_threads - 1) ? ref.size() : (i + 1) * chunk_size;
        threads.emplace_back(buildKmerIndex, ref, std::ref(kmerMap), K, start, end, std::ref(mtx));
    }
    
    for (auto &t : threads) {
        t.join();
    }
    
    // 查找正向匹配的锚点
    for (int i = 0; i + K <= query.size(); ++i) {
        string kmer = query.substr(i, K);
        if (kmerMap.count(kmer)) {
            for (int r_pos : kmerMap[kmer]) {
                anchors.emplace_back(i, i + K, r_pos, r_pos + K);
            }
        }
    }
    
    // 查找反向互补匹配的锚点
    string rc_query = getRC(query);
    for (int i = 0; i + K <= rc_query.size(); ++i) {
        string kmer = rc_query.substr(i, K);
        if (kmerMap.count(kmer)) {
            for (int r_pos : kmerMap[kmer]) {
                int q_start = query.size() - (i + K);
                int q_end = query.size() - i;
                anchors.emplace_back(q_start, q_end, r_pos, r_pos + K, true);
            }
        }
    }
    
    return anchors;
}

vector<tuple<int, int, int, int>> processGaps(const vector<Anchor> &chain, 
                                             const string &query, const string &ref) {
    vector<tuple<int, int, int, int>> result;
    
    if (chain.empty()) {
        result.emplace_back(0, query.size(), 0, ref.size());
        return result;
    }
    
    // 处理第一个锚点之前的区域
    if (chain[0].q_start > 0 || chain[0].r_start > 0) {
        result.emplace_back(0, chain[0].q_start, 0, chain[0].r_start);
    }
    
    // 添加第一个锚点
    result.emplace_back(chain[0].q_start, chain[0].q_end, 
                       chain[0].r_start, chain[0].r_end);
    
    // 处理中间区域和锚点
    for (size_t i = 1; i < chain.size(); ++i) {
        int q_gap_start = chain[i-1].q_end;
        int q_gap_end = chain[i].q_start;
        int r_gap_start = chain[i-1].r_end;
        int r_gap_end = chain[i].r_start;
        
        if (q_gap_start < q_gap_end || r_gap_start < r_gap_end) {
            result.emplace_back(q_gap_start, q_gap_end, r_gap_start, r_gap_end);
        }
        
        result.emplace_back(chain[i].q_start, chain[i].q_end, 
                          chain[i].r_start, chain[i].r_end);
    }
    
    // 处理最后一个锚点之后的区域
    int last_q_end = chain.back().q_end;
    int last_r_end = chain.back().r_end;
    if (last_q_end < query.size() || last_r_end < ref.size()) {
        result.emplace_back(last_q_end, query.size(), last_r_end, ref.size());
    }
    
    return result;
}

string formatOutput(const vector<tuple<int, int, int, int>> &alignment) {
    stringstream ss;
    ss << fixed << setprecision(2);
    ss << "[";
    for (size_t i = 0; i < alignment.size(); ++i) {
        const auto &region = alignment[i];
        ss << "(" << get<0>(region) << "," << get<1>(region) << ","
           << get<2>(region) << "," << get<3>(region) << ")";
        if (i != alignment.size() - 1) {
            ss << ",";
            if (i % 3 == 2) ss << "\n ";
            else ss << " ";
        }
    }
    ss << "]";
    return ss.str();
}

void extendAnchors(vector<Anchor> &anchors, const string &query, const string &ref) {
    mutex mtx;
    
    auto extend = [&](int start, int end) {
        for (int i = start; i < end; ++i) {
            auto &anchor = anchors[i];
            int mismatch_left = 0;
            int ext_left = 0;
            
            // 向左扩展
            while (ext_left < MAX_EXTEND && 
                   anchor.q_start > 0 && 
                   ((!anchor.is_rc && anchor.r_start > 0) || 
                    (anchor.is_rc && anchor.r_end < ref.size()))) {
                char q_base = query[anchor.q_start - 1];
                char r_base = anchor.is_rc ? getRC(string(1, ref[anchor.r_end]))[0] : ref[anchor.r_start - 1];
                
                if (q_base != r_base && q_base != 'N' && r_base != 'N') {
                    if (++mismatch_left > MISMATCH_TOLERANCE) break;
                }
                
                anchor.q_start--;
                if (!anchor.is_rc) {
                    anchor.r_start--;
                } else {
                    anchor.r_end++;
                }
                ext_left++;
            }
            
            int mismatch_right = 0;
            int ext_right = 0;
            // 向右扩展
            while (ext_right < MAX_EXTEND && 
                   anchor.q_end < query.size() && 
                   ((!anchor.is_rc && anchor.r_end < ref.size()) || 
                    (anchor.is_rc && anchor.r_start > 0))) {
                char q_base = query[anchor.q_end];
                char r_base = anchor.is_rc ? getRC(string(1, ref[anchor.r_start - 1]))[0] : ref[anchor.r_end];
                
                if (q_base != r_base && q_base != 'N' && r_base != 'N') {
                    if (++mismatch_right > MISMATCH_TOLERANCE) break;
                }
                
                anchor.q_end++;
                if (!anchor.is_rc) {
                    anchor.r_end++;
                } else {
                    anchor.r_start--;
                }
                ext_right++;
            }
            
            // 重新计算分数
            anchor.score = Anchor::calculateScore(anchor.q_start, anchor.q_end, 
                                                 anchor.r_start, anchor.r_end);
        }
    };
    
    // 并行扩展锚点
    const int num_threads = thread::hardware_concurrency();
    vector<thread> threads;
    int chunk_size = anchors.size() / num_threads;
    
    for (int i = 0; i < num_threads; ++i) {
        int start = i * chunk_size;
        int end = (i == num_threads - 1) ? anchors.size() : (i + 1) * chunk_size;
        threads.emplace_back(extend, start, end);
    }
    
    for (auto &t : threads) {
        t.join();
    }
    
    // 过滤低质量锚点
    anchors.erase(remove_if(anchors.begin(), anchors.end(), 
                    [&](const Anchor &a) { 
                        return (a.q_end - a.q_start) < MIN_ANCHOR_LENGTH || 
                               a.similarity(query, ref) < SIMILARITY_THRESHOLD; 
                    }),
                anchors.end());
}

// 查找未覆盖区域的二级锚点
vector<Anchor> findSecondaryAnchors(const vector<Anchor> &primaryChain, 
                                  const string &query, const string &ref) {
    vector<bool> covered(query.size(), false);
    for (const auto &anchor : primaryChain) {
        for (int i = max(0, anchor.q_start); i < min((int)query.size(), anchor.q_end); ++i) {
            covered[i] = true;
        }
    }
    
    vector<pair<int, int>> uncovered;
    int start = -1;
    for (int i = 0; i < covered.size(); ++i) {
        if (!covered[i]) {
            if (start == -1) start = i;
        } else if (start != -1) {
            uncovered.emplace_back(start, i);
            start = -1;
        }
    }
    if (start != -1) {
        uncovered.emplace_back(start, covered.size());
    }
    
    vector<Anchor> secondaryAnchors;
    int secondaryK = max(K - 4, 6); // 进一步减小k-mer以提高灵敏度
    
    unordered_map<string, vector<int>> kmerMap;
    mutex mtx;
    const int num_threads = thread::hardware_concurrency();
    vector<thread> threads;
    int chunk_size = ref.size() / num_threads;
    
    for (int i = 0; i < num_threads; ++i) {
        int start = i * chunk_size;
        int end = (i == num_threads - 1) ? ref.size() : (i + 1) * chunk_size;
        threads.emplace_back(buildKmerIndex, ref, std::ref(kmerMap), secondaryK, start, end, std::ref(mtx));
    }
    
    for (auto &t : threads) {
        t.join();
    }
    
    for (const auto &region : uncovered) {
        int region_start = region.first;
        int region_end = region.second;
        if (region_end - region_start < secondaryK) continue;
        
        for (int i = region_start; i + secondaryK <= region_end; ++i) {
            string kmer = query.substr(i, secondaryK);
            if (kmerMap.count(kmer)) {
                for (int r_pos : kmerMap[kmer]) {
                    secondaryAnchors.emplace_back(i, i + secondaryK, r_pos, r_pos + secondaryK);
                }
            }
        }
        
        // 反向互补查找
        string rc_query = getRC(query);
        for (int i = region_start; i + secondaryK <= region_end; ++i) {
            string kmer = rc_query.substr(i, secondaryK);
            if (kmerMap.count(kmer)) {
                for (int r_pos : kmerMap[kmer]) {
                    int q_start = query.size() - (i + secondaryK);
                    int q_end = query.size() - i;
                    secondaryAnchors.emplace_back(q_start, q_end, r_pos, r_pos + secondaryK, true);
                }
            }
        }
    }
    
    return secondaryAnchors;
}

int calculateGapPenalty(const Anchor &a1, const Anchor &a2) {
    int q_gap = a2.q_start - a1.q_end;
    int r_gap = a2.is_rc ? (a1.r_start - a2.r_end) : (a2.r_start - a1.r_end);
    int gap_diff = abs(q_gap - r_gap);
    
    // 改进间隔惩罚：基础惩罚 + 非线性缩放
    return GAP_PENALTY_BASE + static_cast<int>(gap_diff * GAP_PENALTY_SCALE * (1 + log2(1 + gap_diff)));
}

vector<Anchor> findOptimalChain(vector<Anchor> &anchors) {
    if (anchors.empty()) return {};
    
    sort(anchors.begin(), anchors.end());
    
    vector<int> dp(anchors.size(), 0);
    vector<int> prev(anchors.size(), -1);
    
    for (int i = 0; i < anchors.size(); ++i) {
        dp[i] = anchors[i].score;
    }
    
    for (int i = 0; i < anchors.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            if (!anchors[j].overlaps(anchors[i]) &&
                anchors[j].q_end <= anchors[i].q_start &&
                ((!anchors[j].is_rc && !anchors[i].is_rc && anchors[j].r_end <= anchors[i].r_start) ||
                 (anchors[j].is_rc && anchors[i].is_rc && anchors[j].r_start >= anchors[i].r_end))) {
                int gap_penalty = calculateGapPenalty(anchors[j], anchors[i]);
                int new_score = dp[j] + anchors[i].score - gap_penalty;
                if (new_score > dp[i]) {
                    dp[i] = new_score;
                    prev[i] = j;
                }
            }
        }
    }
    
    int max_idx = max_element(dp.begin(), dp.end()) - dp.begin();
    vector<Anchor> chain;
    while (max_idx != -1) {
        chain.push_back(anchors[max_idx]);
        max_idx = prev[max_idx];
    }
    
    reverse(chain.begin(), chain.end());
    return chain;
}

vector<tuple<int, int, int, int>> alignSequences(const string &query, const string &ref) {
    // 主锚点查找
    auto anchors = findAnchors(query, ref);
    extendAnchors(anchors, query, ref);
    auto primaryChain = findOptimalChain(anchors);
    
    // 二级锚点查找
    auto secondaryAnchors = findSecondaryAnchors(primaryChain, query, ref);
    if (!secondaryAnchors.empty()) {
        extendAnchors(secondaryAnchors, query, ref);
        primaryChain.insert(primaryChain.end(), secondaryAnchors.begin(), secondaryAnchors.end());
        sort(primaryChain.begin(), primaryChain.end());
        primaryChain = findOptimalChain(primaryChain);
    }
    
    return processGaps(primaryChain, query, ref);
}

int main() {
    string ref, query;
    
    // 从文件读取reference序列
    ifstream ref_file("reference_2_1.txt");
    if (!ref_file) {
        cerr << "无法打开reference.txt文件" << endl;
        return 1;
    }
    getline(ref_file, ref);
    ref_file.close();
    
    // 从文件读取query序列
    ifstream query_file("query_2_1.txt");
    if (!query_file) {
        cerr << "无法打开query.txt文件" << endl;
        return 1;
    }
    getline(query_file, query);
    query_file.close();

    
    auto alignment = alignSequences(query, ref);
    cout << formatOutput(alignment) << endl;

    return 0;
}