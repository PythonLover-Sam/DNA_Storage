# 低级api
def get_reverse_base(base) -> str:
    """
    获得单个互补碱基
    :param base: 需要互补的单个碱基
    :return: 对应的互补碱基
    """
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'C':
        return 'G'
    else:
        return 'C'

def float_to_base(n_float, GC_bias=0.5) -> str:
    """
    给定浮点数转换为碱基 【0，0.25】A 【0.25，0.5】T……
    :param n_float: 浮点数
    :param GC_bias: GC含量偏好
    :return: 碱基
    """
    if n_float < 0:
        raise ValueError("必须大于等于0")
    elif n_float < (1-GC_bias)/2:
        return 'A'
    elif n_float < (1-GC_bias):
        return 'T'
    elif n_float < 1-GC_bias/2:
        return 'C'
    else:
        return 'G'

def get_GC_content(seq) -> float:
    """
    返回一个序列的GC含量
    :param seq: 输入序列
    :return: GC含量
    """
    GC_num = 0
    seq = list(seq)
    total_num = len(seq)
    for base in seq:
        if base == 'G' or base == "C":
            GC_num += 1

    return GC_num / total_num

def simple_get_Tm(seq) -> int:
    """
    计算引物Tm值
    :param seq: 引物序列 长度<=25
    :return:Tm值
    """
    if len(seq) > 25:
        raise ValueError("序列长度必须小于等于25")
    else:
        return int(4 * get_GC_content(seq) * len(seq)
                   + 2 * (1-get_GC_content(seq)) * len(seq))

def check_if_base_complement(base1, base2) -> bool:
    """
    检查两个碱基是否互补
    :param base1: 碱基1
    :param base2: 碱基2
    :return: 是否互补
    """
    if base1 == 'C' and base2 == 'G': return True
    elif base1 == 'G' and base2 == 'C':return True
    elif base1 == 'A' and base2 == 'T': return True
    elif base1 == 'T' and base2 == 'A': return True
    else: return False

def check_max_continues_length(input_list) -> int:
    """
    检查数组中最长的连续数字长度
    :param input_list: 输入数组
    :return: 最大连续长度
    """
    max_ctn_cmpl_num = 0
    if len(input_list) > 0:
        _max_num = 1
        tmp = input_list[0]
        for x in range(len(input_list) - 1):
            if tmp + 1 == input_list[x + 1]:
                tmp += 1
                _max_num += 1
            else:
                tmp = input_list[x + 1]
                _max_num = 1
            if max_ctn_cmpl_num < _max_num:
                max_ctn_cmpl_num = _max_num

        return max_ctn_cmpl_num
    else:
        return 0

def check_max_end_continues_length(input_list, primer_length=20, consider_pos_span=0) -> int:
    """
    检查序列末端在给定跨度下的最大连续互补碱基数
    :param input_list: 输入序列
    :param consider_pos_span: 考虑的跨度 0为严格从末尾最后一个开始向前计算 1为考虑末尾最后一个和倒数第二个开始 以此类推
    :return: 最大末尾连续互补碱基数
    """
    l = list(input_list)
    max_end_continues_length = 0
    m = primer_length - 1
    for i in range(consider_pos_span+1):
        end_continues_length = 1
        m -= 1
        for j in range(len(l) - i - 1):
            if m == l[len(l) - i - 2 - j]:
                m -= 1
                end_continues_length += 1
            else:

                break
        if max_end_continues_length < end_continues_length:
            max_end_continues_length = end_continues_length

    return max_end_continues_length

def check_hamming_distance_between_two_seqs(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("两序列长度不一致")
    else:
        distance = 0
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]:
                distance += 1
        return distance