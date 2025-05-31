# DNA-PJ-Lab2

## 1. 算法概述

本实验旨在分析给定的 DNA 参考序列和查询序列，识别其中的重复片段，并对其进行相似性分析。基于 k-mer锚点 和 滑动窗口 的比对算法。。

## 2. 算法伪代码

### 2.1 主比对算法

```
ALGORITHM DNASequenceAlignment(ref, query)
INPUT: 参考序列ref，查询序列query
OUTPUT: 比对结果列表alignments

BEGIN
    all_alignments ← []
    
    // 第一轮：基于k-mer的比对
    FOR each (kmer_size, step) in kmer_params DO
        anchors ← FindAnchors(ref, query, kmer_size, step)
        IF anchors is not empty THEN
            chain ← ChainAnchors(anchors, ref, query)
            refined_chain ← RefineAlignment(ref, query, chain)
            all_alignments.extend(refined_chain)
        END IF
    END FOR
    
    // 第二轮：滑动窗口比对
    window_params ← SelectWindowSize(query)
    window_alignments ← SlidingWindowAlign(ref, query, window_params)
    all_alignments.extend(window_alignments)
    
    // 合并和去重
    RETURN MergeAlignments(all_alignments)
END
```

### 2.2 k-mer锚点查找算法

```
ALGORITHM FindAnchors(ref, query, kmer_size, step)
INPUT: 参考序列ref，查询序列query，k-mer大小kmer_size，步长step
OUTPUT: 锚点列表anchors

BEGIN
    ref_kmers ← BuildKmerIndex(ref, kmer_size)
    query_rc ← GetReverseComplement(query)
    anchors ← []
    
    FOR i ← 0 TO len(query) - kmer_size STEP step DO
        // 正向匹配
        kmer ← query[i:i+kmer_size]
        IF kmer in ref_kmers THEN
            FOR each ref_pos in ref_kmers[kmer] DO
                anchors.append((i, ref_pos, 1, kmer_size))
            END FOR
        END IF
        
        // 反向互补匹配
        rc_kmer ← query_rc[len(query_rc)-(i+kmer_size):len(query_rc)-i]
        IF rc_kmer in ref_kmers THEN
            FOR each ref_pos in ref_kmers[rc_kmer] DO
                anchors.append((i, ref_pos, -1, kmer_size))
            END FOR
        END IF
    END FOR
    
    RETURN anchors
END
```

### 2.3 锚点扩展算法

```
ALGORITHM ExtendAnchor(ref, query, anchor)
INPUT: 参考序列ref，查询序列query，锚点anchor
OUTPUT: 扩展后的锚点extended_anchor

BEGIN
    (query_pos, ref_pos, strand, length) ← anchor
    
    IF strand == 1 THEN  // 正向扩展
        (left_ext, right_ext) ← ExtendForward(ref, query, query_pos, ref_pos, length)
        RETURN (query_pos - left_ext, query_pos + length + right_ext,
                ref_pos - left_ext, ref_pos + length + right_ext)
    ELSE  // 反向互补扩展
        query_rc ← GetReverseComplement(query)
        rc_query_pos ← len(query) - query_pos - length
        (left_ext, right_ext) ← ExtendForward(ref, query_rc, rc_query_pos, ref_pos, length)
        RETURN (query_pos - right_ext, query_pos + length + left_ext,
                ref_pos - left_ext, ref_pos + length + right_ext)
    END IF
END
```

### 2.4 动态规划链接锚点算法

```
ALGORITHM ChainAnchorsDP(anchors)
INPUT: 锚点列表anchors
OUTPUT: 最优锚点链chain

BEGIN
    n ← len(anchors)
    dp ← array[n] initialized to 0
    parent ← array[n] initialized to -1
    
    FOR i ← 0 TO n-1 DO
        (q_start, q_end, r_start, r_end) ← anchors[i]
        dp[i] ← q_end - q_start  // 初始分数为匹配长度
        
        FOR j ← 0 TO i-1 DO
            (prev_q_start, prev_q_end, prev_r_start, prev_r_end) ← anchors[j]
            
            // 检查是否可以连接
            IF prev_q_end <= q_start AND CanConnect(prev_anchor, current_anchor) THEN
                score ← dp[j] + (q_end - q_start)
                IF score > dp[i] THEN
                    dp[i] ← score
                    parent[i] ← j
                END IF
            END IF
        END FOR
    END FOR
    
    // 回溯找最优路径
    best_idx ← argmax(dp)
    path ← []
    idx ← best_idx
    WHILE idx != -1 DO
        path.append(anchors[idx])
        idx ← parent[idx]
    END WHILE
    
    RETURN reverse(path)
END
```

### 2.5 滑动窗口比对算法

```
ALGORITHM SlidingWindowAlign(ref, query, window_size, step)
INPUT: 参考序列ref，查询序列query，窗口大小window_size，步长step
OUTPUT: 比对结果列表alignments

BEGIN
    alignments ← []
    is_long ← IsLongSequence(query)
    step_size ← SelectStepSize(is_long)
    threshold ← SelectThreshold(is_long)
    
    FOR q_start ← 0 TO len(query) - window_size STEP step DO
        q_end ← min(q_start + window_size, len(query))
        query_window ← query[q_start:q_end]
        query_rc ← GetReverseComplement(query_window)
        
        best_score ← -1
        best_match ← None
        
        FOR r_start ← 0 TO len(ref) - len(query_window) STEP step_size DO
            r_end ← r_start + len(query_window)
            ref_window ← ref[r_start:r_end]
            
            // 计算正向和反向互补的编辑距离
            dist1 ← EditDistance(ref_window, query_window)
            dist2 ← EditDistance(ref_window, query_rc)
            
            score1 ← len(query_window) - dist1
            score2 ← len(query_window) - dist2
            score ← max(score1, score2)
            
            IF score > best_score AND score >= len(query_window) * threshold THEN
                best_score ← score
                best_match ← (q_start, q_end, r_start, r_end)
            END IF
        END FOR
        
        IF best_match is not None THEN
            alignments.append(best_match)
        END IF
    END FOR
    
    RETURN alignments
END
```

## 3. 时空复杂度分析

### 3.1 时间复杂度

设参考序列长度为n，查询序列长度为m。
- 使用动态规划和k-mer哈希方法
对于实际参数：
- 长序列：w=498, s=1，复杂度约为O(mn)
- 短序列：w=100, s=4，复杂度约为O(mn/4)

### 3.2 空间复杂度

- **k-mer索引**：O(n)，存储所有k-mer及其位置
- **锚点存储**：O(k)，存储所有找到的锚点
- **动态规划数组**：O(k)，存储DP状态
- **比对结果**：O(m/w)，存储滑动窗口比对结果

**总空间复杂度**：O(n + k + m/w) = O(n + m)


## 4. 运行结果

### 4.1 测试数据
- **输入文件**：text1.txt和text2.txt
- **数据格式**：包含参考序列(ref)和查询序列(query)

### 4.2 比对结果

程序成功运行，输出了大量比对片段，结果格式为：
```
[(query_start, query_end, ref_start, ref_end), ...]
```

## 5. 

短的处理成功了，结果2107，长的还是差一点（29.806），实在调不出来，放弃了。