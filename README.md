# dna-pj
# DNA 序列重复片段查找实验报告

## 1. 实验目的

本实验旨在分析给定的 DNA 参考序列和查询序列，识别其中的重复片段，并对其进行相似性分析。同时，考虑 DNA 互补链的重复情况。

## 2. 实验方法

本实验使用动态规划方法计算最长公共子序列（LCS），以识别 DNA 片段的重复情况。同时，通过合并相似的重复片段提高分析的准确性。

## 3. 伪代码

### 读取 DNA 序列

```python
Function read_file(file_path):
    Open file at file_path
    Read and return stripped content
```

### 计算互补序列

```python
Function get_complement_sequence(sequence):
    Define complement dictionary {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    Reverse traverse sequence:
        Append complement base to complement_seq
    Return complement_seq as string
```

### 计算最长公共子序列（LCS）

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

### 查找重复片段

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

### 合并相似重复片段

```python
Function merge_similar_repeats(repeats, length_diff=2, similarity_threshold=0.95):
    Group repeats by complement status
    Sort groups by segment length
    Merge similar repeats if:
        Length difference <= length_diff
        Similarity >= similarity_threshold
    Return merged repeats
```

## 4. 复杂度分析

### 读取文件

**时间复杂度**：O(N)，N 为 DNA 序列长度。 **空间复杂度**：O(1)。

### 计算互补序列

**时间复杂度**：O(N) **空间复杂度**：O(N)

### 计算最长公共子序列（LCS）

**时间复杂度**：O(M × N) **空间复杂度**：O(M × N) 其中 M 和 N 分别为参考序列和查询序列的长度。

### 查找重复片段

**时间复杂度**：O(M × N) **空间复杂度**：O(M × N)

### 合并相似重复片段

**时间复杂度**：O(K^2)，K 为重复片段数量（假设使用冒泡排序）。 **空间复杂度**：O(K)

## 5. 实验总结

本实验成功实现了 DNA 序列重复片段的识别和合并，通过 LCS 方法进行相似性计算，并综合考虑了反向互补序列的重复情况。未来可优化排序算法和内存占用，提高处理大规模数据的效率。
