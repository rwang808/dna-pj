import numpy as np
import edlib
from collections import defaultdict
from typing import List, Tuple, Optional
import heapq


class DNASequenceAligner:
    """DNA序列比对器类，封装所有比对相关的功能"""
    
    def __init__(self):
        """初始化比对器配置"""
        # K-mer参数配置
        self.kmer_params = [
            (12, 3),
            (15, 5),
            (10, 2),
            (8, 1),
            (14, 4),
            (16, 6),
            (9, 1),
        ]
        
        # 序列长度阈值
        self.long_sequence_threshold = 10000
        
        # 扩展参数
        self.max_extension = 50
        self.min_anchor_length = 15
        
        # 精化参数
        self.long_seq_edit_threshold = 0.05
        self.short_seq_edit_threshold = 0.16
        self.min_match_length = 30
        self.split_threshold = 50
        
        # 滑动窗口参数
        self.long_seq_window = (498, 498)
        self.short_seq_window = (100, 100)
        self.long_seq_step_size = 1
        self.short_seq_step_size = 4
        self.long_seq_threshold = 0.5
        self.short_seq_threshold = 0.52
        
        # 链接参数
        self.max_gap_distance = 2000
        
        # 互补碱基映射
        self.complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    
    def get_reverse_complement(self, sequence: str) -> str:
        """获取反向互补序列"""
        complement = [self.complement_map[base] for base in sequence]
        return ''.join(complement[::-1])
    
    def is_long_sequence(self, query: str) -> bool:
        """判断是否为长序列"""
        return len(query) > self.long_sequence_threshold
    
    def build_kmer_index(self, ref: str, kmer_size: int) -> defaultdict:
        """构建k-mer索引"""
        ref_kmers = defaultdict(list)
        for i in range(len(ref) - kmer_size + 1):
            kmer = ref[i:i + kmer_size]
            if 'N' not in kmer:
                ref_kmers[kmer].append(i)
        return ref_kmers
    
    def find_anchors(self, ref: str, query: str, kmer_size: int = 12, step: int = 3) -> List[Tuple]:
        """使用k-mer哈希找到锚点"""
        ref_kmers = self.build_kmer_index(ref, kmer_size)
        query_rc = self.get_reverse_complement(query)
        anchors = []
        
        # 在查询序列中寻找匹配的k-mer
        for i in range(0, len(query) - kmer_size + 1, step):
            # 正向匹配
            kmer = query[i:i + kmer_size]
            if kmer in ref_kmers:
                for ref_pos in ref_kmers[kmer]:
                    anchors.append((i, ref_pos, 1, kmer_size))
            
            # 反向互补匹配
            rc_kmer = query_rc[len(query_rc) - (i + kmer_size):len(query_rc) - i]
            if rc_kmer in ref_kmers:
                for ref_pos in ref_kmers[rc_kmer]:
                    anchors.append((i, ref_pos, -1, kmer_size))
        
        return anchors
    
    def extend_anchor(self, ref: str, query: str, anchor: Tuple) -> Tuple:
        """扩展锚点以获得更长的匹配"""
        query_pos, ref_pos, strand, length = anchor
        
        if strand == 1:
            # 正向扩展
            left_ext, right_ext = self._extend_forward(ref, query, query_pos, ref_pos, length)
            return (query_pos - left_ext, query_pos + length + right_ext,
                   ref_pos - left_ext, ref_pos + length + right_ext)
        else:
            # 反向互补扩展
            query_rc = self.get_reverse_complement(query)
            rc_query_pos = len(query) - query_pos - length
            left_ext, right_ext = self._extend_forward(ref, query_rc, rc_query_pos, ref_pos, length)
            return (query_pos - right_ext, query_pos + length + left_ext,
                   ref_pos - left_ext, ref_pos + length + right_ext)
    
    def _extend_forward(self, ref: str, query: str, query_pos: int, ref_pos: int, length: int) -> Tuple[int, int]:
        """正向扩展辅助函数"""
        left_ext = 0
        right_ext = 0
        
        # 向左扩展
        while (query_pos - left_ext - 1 >= 0 and
               ref_pos - left_ext - 1 >= 0 and
               left_ext < self.max_extension and
               query[query_pos - left_ext - 1] == ref[ref_pos - left_ext - 1]):
            left_ext += 1
        
        # 向右扩展
        while (query_pos + length + right_ext < len(query) and
               ref_pos + length + right_ext < len(ref) and
               right_ext < self.max_extension and
               query[query_pos + length + right_ext] == ref[ref_pos + length + right_ext]):
            right_ext += 1
        
        return left_ext, right_ext
    
    def calculate_edit_distance(self, ref: str, query: str, ref_start: int, ref_end: int, 
                              query_start: int, query_end: int) -> int:
        """计算两个序列片段之间的编辑距离"""
        ref_seq = ref[ref_start:ref_end]
        query_seq = query[query_start:query_end]
        query_rc = self.get_reverse_complement(query_seq)
        
        dist1 = edlib.align(ref_seq, query_seq)['editDistance']
        dist2 = edlib.align(ref_seq, query_rc)['editDistance']
        
        return min(dist1, dist2)
    
    def merge_overlapping_anchors(self, anchors: List[Tuple]) -> List[Tuple]:
        """合并重叠的锚点"""
        if not anchors:
            return []
        
        anchors.sort(key=lambda x: (x[0], x[1]))
        merged_anchors = []
        
        for anchor in anchors:
            q_start, q_end, r_start, r_end = anchor
            merged = False
            
            for i, (mq_start, mq_end, mr_start, mr_end) in enumerate(merged_anchors):
                # 检查是否重叠
                if (q_start < mq_end and q_end > mq_start and
                        r_start < mr_end and r_end > mr_start):
                    # 合并重叠区域，选择更长的
                    if (q_end - q_start) > (mq_end - mq_start):
                        merged_anchors[i] = anchor
                    merged = True
                    break
            
            if not merged:
                merged_anchors.append(anchor)
        
        return sorted(merged_anchors, key=lambda x: x[0])
    
    def chain_anchors_dp(self, anchors: List[Tuple]) -> List[Tuple]:
        """使用动态规划链接锚点"""
        if not anchors:
            return []
        
        n = len(anchors)
        dp = [0] * n
        parent = [-1] * n
        
        # 动态规划找最优链
        for i in range(n):
            q_start, q_end, r_start, r_end = anchors[i]
            dp[i] = q_end - q_start  # 初始分数为匹配长度
            
            for j in range(i):
                prev_q_start, prev_q_end, prev_r_start, prev_r_end = anchors[j]
                
                # 检查是否可以连接（无重叠且顺序合理）
                if (prev_q_end <= q_start and
                        (prev_r_end <= r_start or 
                         abs((q_start - prev_q_end) - (r_start - prev_r_end)) < self.max_gap_distance)):
                    
                    score = dp[j] + (q_end - q_start)
                    if score > dp[i]:
                        dp[i] = score
                        parent[i] = j
        
        # 回溯找最优路径
        best_idx = np.argmax(dp)
        path = []
        idx = best_idx
        while idx != -1:
            path.append(anchors[idx])
            idx = parent[idx]
        
        path.reverse()
        return path
    
    def chain_anchors(self, anchors: List[Tuple], ref: str, query: str) -> List[Tuple]:
        """链接锚点的主函数"""
        if not anchors:
            return []
        
        # 扩展锚点
        extended_anchors = []
        for anchor in anchors:
            try:
                ext = self.extend_anchor(ref, query, anchor)
                if ext[1] - ext[0] >= self.min_anchor_length:
                    extended_anchors.append(ext)
            except Exception:
                continue
        
        if not extended_anchors:
            return []
        
        # 合并重叠的锚点
        merged_anchors = self.merge_overlapping_anchors(extended_anchors)
        
        # 使用动态规划链接
        return self.chain_anchors_dp(merged_anchors)
    
    def split_segment(self, ref: str, query: str, segment: Tuple) -> List[Tuple]:
        """分割长片段"""
        q_start, q_end, r_start, r_end = segment
        mid_q = (q_start + q_end) // 2
        mid_r = (r_start + r_end) // 2
        
        parts = [
            (q_start, mid_q, r_start, mid_r),
            (mid_q, q_end, mid_r, r_end)
        ]
        
        refined_parts = []
        edit_threshold = (self.long_seq_edit_threshold if self.is_long_sequence(query) 
                         else self.short_seq_edit_threshold)
        
        for part in parts:
            pq_start, pq_end, pr_start, pr_end = part
            if pq_end - pq_start >= self.min_match_length:
                pedit_dist = self.calculate_edit_distance(ref, query, pr_start, pr_end, pq_start, pq_end)
                if pedit_dist / (pq_end - pq_start) <= edit_threshold:
                    refined_parts.append(part)
        
        return refined_parts
    
    def refine_alignment(self, ref: str, query: str, initial_chain: List[Tuple]) -> List[Tuple]:
        """精化比对结果"""
        if not initial_chain:
            return []
        
        refined = []
        edit_threshold = (self.long_seq_edit_threshold if self.is_long_sequence(query) 
                         else self.short_seq_edit_threshold)
        
        for segment in initial_chain:
            q_start, q_end, r_start, r_end = segment
            edit_dist = self.calculate_edit_distance(ref, query, r_start, r_end, q_start, q_end)
            match_len = q_end - q_start
            
            if match_len >= self.min_match_length and edit_dist / match_len <= edit_threshold:
                refined.append(segment)
            elif match_len >= self.split_threshold and edit_dist / match_len > edit_threshold:
                # 尝试分割成更小的片段
                refined.extend(self.split_segment(ref, query, segment))
        
        return refined
    
    def sliding_window_align(self, ref: str, query: str, window_size: int = 50, step: int = 25) -> List[Tuple]:
        """滑动窗口比对，用于填补空隙"""
        alignments = []
        is_long = self.is_long_sequence(query)
        step_size = self.long_seq_step_size if is_long else self.short_seq_step_size
        threshold = self.long_seq_threshold if is_long else self.short_seq_threshold
        
        for q_start in range(0, len(query) - window_size + 1, step):
            q_end = min(q_start + window_size, len(query))
            query_window = query[q_start:q_end]
            query_rc = self.get_reverse_complement(query_window)
            
            best_score = -1
            best_match = None
            
            # 在参考序列中搜索最佳匹配
            for r_start in range(0, len(ref) - len(query_window) + 1, step_size):
                r_end = r_start + len(query_window)
                ref_window = ref[r_start:r_end]
                
                # 正向和反向互补匹配
                dist1 = edlib.align(ref_window, query_window)['editDistance']
                dist2 = edlib.align(ref_window, query_rc)['editDistance']
                
                score1 = len(query_window) - dist1
                score2 = len(query_window) - dist2
                score = max(score1, score2)
                
                if score > best_score and score >= len(query_window) * threshold:
                    best_score = score
                    best_match = (q_start, q_end, r_start, r_end)
            
            if best_match:
                alignments.append(best_match)
        
        return alignments
    
    def merge_alignments(self, alignments: List[Tuple]) -> List[Tuple]:
        """合并和去重比对结果"""
        if not alignments:
            return []
        
        # 按query位置排序
        alignments.sort(key=lambda x: x[0])
        
        # 去重并合并重叠区域
        merged = []
        for alignment in alignments:
            q_start, q_end, r_start, r_end = alignment
            merged_flag = False
            
            for i, (mq_start, mq_end, mr_start, mr_end) in enumerate(merged):
                # 检查重叠
                if q_start < mq_end and q_end > mq_start:
                    # 选择更长的片段
                    if (q_end - q_start) > (mq_end - mq_start):
                        merged[i] = alignment
                    merged_flag = True
                    break
            
            if not merged_flag:
                merged.append(alignment)
        
        return merged
    
    def align_sequences(self, ref: str, query: str) -> List[Tuple]:
        """主要的序列比对函数"""
        all_alignments = []
        
        # 第一轮：基于k-mer的比对
        for kmer_size, step in self.kmer_params:
            anchors = self.find_anchors(ref, query, kmer_size=kmer_size, step=step)
            
            if not anchors:
                continue
            
            # 链接锚点
            chain = self.chain_anchors(anchors, ref, query)
            
            # 精化比对
            refined_chain = self.refine_alignment(ref, query, chain)
            all_alignments.extend(refined_chain)
        
        # 第二轮：根据序列长度自适应选择窗口尺寸
        window_params = (self.long_seq_window if self.is_long_sequence(query) 
                        else self.short_seq_window)
        
        window_alignments = self.sliding_window_align(ref, query, 
                                                     window_size=window_params[0], 
                                                     step=window_params[1])
        all_alignments.extend(window_alignments)
        
        # 合并和去重
        return self.merge_alignments(all_alignments)
    
    def calculate_alignment_score(self, ref: str, query: str, alignment: Tuple) -> int:
        """计算比对片段的得分"""
        q_start, q_end, r_start, r_end = alignment
        match_len = q_end - q_start
        edit_dist = self.calculate_edit_distance(ref, query, r_start, r_end, q_start, q_end)
        return max(0, match_len - edit_dist)
    
    def post_process_alignments(self, ref: str, query: str, alignments: List[Tuple]) -> List[Tuple]:
        """后处理比对结果，进行高质量去重"""
        unique_alignments = []
        
        for alignment in alignments:
            q_start, q_end, r_start, r_end = alignment
            merged = False
            
            for i, existing in enumerate(unique_alignments):
                eq_start, eq_end, er_start, er_end = existing
                
                # 检查query位置重叠
                if not (q_end <= eq_start or q_start >= eq_end):
                    # 计算当前和已存在片段的质量
                    current_score = self.calculate_alignment_score(ref, query, alignment)
                    existing_score = self.calculate_alignment_score(ref, query, existing)
                    
                    if current_score > existing_score:
                        unique_alignments[i] = alignment
                    merged = True
                    break
            
            if not merged:
                unique_alignments.append(alignment)
        
        return unique_alignments


def format_output(alignments: List[Tuple]) -> str:
    """格式化输出"""
    result = []
    for q_start, q_end, r_start, r_end in alignments:
        result.append(f"({q_start}, {q_end}, {r_start}, {r_end})")
    return "[" + ", ".join(result) + "]"


def load_test_data(file_path: str) -> Tuple[str, str]:
    """加载测试数据"""
    with open(file_path, 'r') as f:
        content = f.read().strip()
    
    # 解析数据
    parts = content.split('\n\nquery:\n')
    ref = parts[0].split('\n')[1]  # 跳过第一行标签
    query = parts[1]
    
    return ref, query


def main():
    """主函数"""
    # 创建比对器实例
    aligner = DNASequenceAligner()

    ref, query = load_test_data(r"D:\lab\choice\text1.txt")
    
    # 执行比对
    alignments = aligner.align_sequences(ref, query)
    
    # 去重
    unique_alignments = aligner.post_process_alignments(ref, query, alignments)

    print(format_output(unique_alignments))



if __name__ == "__main__":
    main()