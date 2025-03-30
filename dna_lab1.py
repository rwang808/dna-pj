
def read_file(file_path):
    with open(file_path, 'r') as f:
        return f.read().strip()


def get_complement_sequence(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement_seq = []
    for i in range(len(sequence) - 1, -1, -1):
        complement_seq.append(complement.get(sequence[i], sequence[i]))
    return ''.join(complement_seq)

def sequence_similarity(seq1, seq2):
    """简易相似度计算"""
    len1 = len(seq1)
    len2 = len(seq2)

    # 计算最长公共子序列长度
    dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if seq1[i - 1] == seq2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])

    lcs_length = dp[len1][len2]
    max_len = max(len1, len2)
    return lcs_length / max_len if max_len > 0 else 0.0


def find_repeats(reference, query, is_complement=False):
    ref_len = len(reference)
    query_len = len(query)
    dp = [[0] * (query_len + 1) for _ in range(ref_len + 1)]
    repeats = []
    threshold = 20

    for i in range(1, ref_len + 1):
        for j in range(1, query_len + 1):
            if reference[i - 1] == query[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
                if dp[i][j] >= threshold and (i == ref_len or j == query_len or reference[i] != query[j]):
                    length = dp[i][j]
                    ref_start = i - length
                    repeats.append({
                        'length': length,
                        'ref_position': ref_start + length,
                        'segment': reference[ref_start:i],
                        'is_complement': is_complement,
                        'repeat_count': 1
                    })
    return repeats


def merge_similar_repeats(repeats, length_diff=2, similarity_threshold=0.95):
    if not repeats:
        return []

    groups = {}
    for repeat in repeats:
        key = (repeat['is_complement'])
        if key not in groups:
            groups[key] = []
        groups[key].append(repeat)

    merged = []
    for group in groups.values():
        # 按长度排序 (冒泡排序)
        n = len(group)
        for i in range(n):
            for j in range(0, n - i - 1):
                if group[j]['length'] > group[j + 1]['length']:
                    group[j], group[j + 1] = group[j + 1], group[j]

        i = 0
        while i < len(group):
            current = group[i]
            j = i + 1
            while j < len(group):
                sim = sequence_similarity(current['segment'], group[j]['segment'])
                if (abs(current['length'] - group[j]['length']) <= length_diff and
                        sim >= similarity_threshold):

                    current['repeat_count'] += group[j]['repeat_count']
                    current['length'] = max(current['length'], group[j]['length'])
                    current['ref_position'] = (current['ref_position'] + group[j]['ref_position']) // 2
                    j += 1
                else:
                    break

            merged.append(current)
            i = j

    return merged




def print_combined_results(normal_repeats, comp_repeats):
    threshold = 20
    all_repeats = normal_repeats + comp_repeats
    filtered = []
    for r in all_repeats:
        if r['length'] >= threshold and r['repeat_count'] > 1:
            filtered.append(r)

    # 排序：长度降序，位置升序 (冒泡排序)
    n = len(filtered)
    for i in range(n):
        for j in range(0, n - i - 1):
            swap = False
            if filtered[j]['length'] < filtered[j + 1]['length']:
                swap = True
            elif (filtered[j]['length'] == filtered[j + 1]['length'] and
                  filtered[j]['ref_position'] > filtered[j + 1]['ref_position']):
                swap = True

            if swap:
                filtered[j], filtered[j + 1] = filtered[j + 1], filtered[j]

    print("+----------+-----------+------------+----------+-----------+----------+")
    print("| 序号     | 重复片段   | pos in ref  | 序列长度  | 重复数量   | 是否反向  |")
    print("+----------+-----------+------------+----------+-----------+----------+")

    for i in range(len(filtered)):
        r = filtered[i]
        seg = r['segment'][:20] + "..." if len(r['segment']) > 20 else r['segment']
        direction = "是" if r['is_complement'] else "否"
        print(f"| {i + 1:<8} | {seg:<3} | {r['ref_position']:<4} | "
              f"{r['length']:<4} | {r['repeat_count']:<4} | "
              f"{direction:<4} |")

    print("+----------+-----------+------------+----------+-----------+----------+")


def main():
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    ref = read_file(os.path.join(current_dir, 'ref.txt'))
    query = read_file(os.path.join(current_dir, 'query.txt'))

    print(f"参考序列长度: {len(ref)}")
    print(f"查询序列长度: {len(query)}")

    # 处理正向重复
    normal = find_repeats(ref, query)
    normal = merge_similar_repeats(normal)

    # 处理反向互补重复
    comp_query = get_complement_sequence(query)
    comp_repeats = find_repeats(ref, comp_query, True)
    comp_repeats = merge_similar_repeats(comp_repeats)

    # 打印结果
    print_combined_results(normal, comp_repeats)


if __name__ == "__main__":
    main()