#include <vector>
#include <string>
#include <tuple>
#include <iostream>
#include <unordered_map>
#include <algorithm>

using namespace std;

// 获取反向互补序列
string getRC(const string &s)
{
    string rc;
    for (char c : s) 
    {
        switch (c) 
        {
            case 'A': rc += 'T'; break;
            case 'T': rc += 'A'; break;
            case 'C': rc += 'G'; break;
            case 'G': rc += 'C'; break;
        }
    }

    reverse(rc.begin(), rc.end());
    return rc;
}


//填充DP表并查找重复片段 ----- 使用路径规划实现
void finding(const string &ref, const string &query, unordered_map<string, bool> &reported,
                    vector<tuple<int, int, int, bool>> &repeats, bool isRC) 
{
    int n = ref.size();
    int m = query.size();
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));

    for (int i = 1; i <= n; i++) 
    {
        for (int j = 1; j <= m; j++) 
        {
            // 填充DP表
            if (ref[i-1] == query[j-1]) 
            {
                dp[i][j] = dp[i-1][j-1] + 1;
                
                int len = dp[i][j];
                if (len >= 10) 
                {
                    string sub = ref.substr(i - len, len);
                    string key = isRC ? getRC(sub) : sub;
                    
                    // 确保不重复报告
                    if (!reported[key]) 
                    {
                        // 计算连续重复次数
                        int repeat_count = 1;
                        int pos = j + 1;  // 从下一个位置开始检查
                        
                        while (pos + len - 1 <= m && query.substr(pos - 1, len) == sub) 
                        {
                            repeat_count++;
                            pos += len;
                        }
                        
                        // 只记录有重复的片段
                        if (repeat_count > 1) 
                        {
                            repeats.emplace_back(i, len, repeat_count, isRC);
                            reported[key] = true;
                        }
                    }
                }
            } 
            else 
            {
                dp[i][j] = 0;
            }
        }
    }
}

// 快速排序分区函数
int partition(vector<tuple<int, int, int, bool>>& arr, int low, int high) {
    int pivot = get<1>(arr[high]);  // 选择最后一个元素作为基准
    int i = low - 1;                // 小于pivot的元素的索引

    for (int j = low; j <= high - 1; j++) 
    {
        // 按长度降序排列
        if (get<1>(arr[j]) < pivot) 
        {
            i++;
            swap(arr[i], arr[j]);
        }
    }
    swap(arr[i + 1], arr[high]);
    return i + 1;
}

// 自定义快速排序
void quickSort(vector<tuple<int, int, int, bool>>& arr, int low, int high) 
{
    if (low < high) {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

// 主函数：寻找重复片段
vector<tuple<int, int, int, bool>> findRepeats(const string &query, const string &ref) 
{
    vector<tuple<int, int, int, bool>> repeats;
    unordered_map<string, bool> reported;
    
    finding(ref, query, reported, repeats, false);
    
    string rc_query = getRC(query);
    finding(ref, rc_query, reported, repeats, true);
    
    //用自定义快排进行排序 --- 结果从小到大输出
    if (!repeats.empty()) 
    {
        quickSort(repeats, 0, repeats.size() - 1);
    }
    
    return repeats;
}


int main() 
{
    string ref, query;

    cout << "input reference:" << endl;
    cin >> ref;

    cout << "input query:" << endl;
    cin >> query;

    auto repeats = findRepeats(query, ref);
    
    // 输出结果
    for (const auto &repeat : repeats) {
        cout << "POS IN REF: " << get<0>(repeat)<< endl;  
        cout << "REPEAT SIZE: " << get<1>(repeat) << endl;
        cout << "REPEAT COUNT: " << get<2>(repeat) << endl;
        cout << "INVERSE: " << (get<3>(repeat) ? "True" : "False") << endl;
        cout << "--------------------" << endl;
    }
    
    return 0;
}