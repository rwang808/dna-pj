# dna-pj
# DNA 序列重复片段lab1

## 1. 实验目的

本实验旨在分析给定的 DNA 参考序列和查询序列，识别其中的重复片段，并对其进行相似性分析。同时，考虑 DNA 互补链的重复情况。

## 2. 实验方法

本实验使用动态规划方法计算最长公共子序列，以识别 DNA 片段的重复情况。同时，通过合并相似的重复片段提高分析的准确性。

## 3. 算法伪代码


### 计算互补序列

```python
Function get_complement_sequence(sequence):
    Define complement dictionary {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    Reverse traverse sequence:
        Append complement base to complement_seq
    Return complement_seq as string
```

### 计算最长公共子序列

```python
Function sequence_similarity(seq1, seq2):
    Initialize dp table of size (len(seq1) + 1) x (len(seq2) + 1)
    For each character in seq1:
        For each character in seq2:
            If characters match:
                dp[i][j] = dp[i-1][j-1] + 1
            Else:
                dp[i][j] = max(dp[i-1][j], dp[i][j-1])
    Return dp[len(seq1)][len(seq2)] / max(len(seq1), len(seq2))
```

### 动态规划查找重复片段

```python
Function find_repeats(reference, query, is_complement=False):
    Initialize dp table of size (len(reference) + 1) x (len(query) + 1)
    Initialize empty repeat list
    Define threshold = 20
    For each character in reference:
        For each character in query:
            If characters match:
                dp[i][j] = dp[i-1][j-1] + 1
                If dp[i][j] >= threshold and end of sequence:
                    Store repeat information
    Return repeat list
```


## 4. 复杂度分析

**时间复杂度**：使用动态规划方法寻找两个序列中的重复片段时间复杂度：O(n²)，其中n是序列长度。 

**空间复杂度**：也是O(n²)。


## 5. 运行结果



