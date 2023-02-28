# DNA序列处理库
import Functions.LowLevelAPI as lapi
import random
import primer3
from copy import deepcopy
from Parameters import Parameter
from progressbar import *


def get_reverse_strand(seq) -> str:
    """
    获取输入序列的反向序列
    :param seq:
    :return:
    """
    l = list(seq)
    l.reverse()
    return ''.join(l)

def get_complement_strand(seq, if_5To3=True) -> str:
    """
    获取给定序列的互补序列，需要输入5‘ - 3’顺序
    :param seq: 输入序列
    :param if_5To3:是否按照5‘ - 3’顺序
    :return: 互补序列
    """
    l_seq = list(seq)
    for i in range(len(l_seq)):
        l_seq[i] = lapi.get_reverse_base(l_seq[i])
    if if_5To3:
        l_seq.reverse()
        return ''.join(l_seq)
    return ''.join(l_seq)

def create_random_strand(length, GC_content=0.5) -> str:
    """
    创建随机序列
    :param length:序列长度
    :param GC_content: 偏好GC含量
    :return: 随机序列
    """
    l_rand_num = []
    for _ in range(length):
        l_rand_num.append(random.random())
    l_rand_base = []
    for item in l_rand_num:
        l_rand_base.append(lapi.float_to_base(item, GC_content))
    return ''.join(l_rand_base)

def check_if_satisfy_primer_principle(seq, consider_tm=True) -> int:
    if seq[-1] == 'A':
        return -5
    if lapi.get_GC_content(seq) < Parameter.primer_GC_content_min or \
            lapi.get_GC_content(seq) > Parameter.primer_GC_content_max:
        return -1
    if consider_tm:
        if abs(lapi.simple_get_Tm(''.join(seq)) - 57.5) >= Parameter.primer_Tm_bias_from_57_5_degree:
            return -2
    if primer3.calcHairpin(''.join(seq)).structure_found is True:
        return -3
    if seq[-2] == seq[-3] == seq[-4] == 'C' or \
            seq[-2] == seq[-3] == seq[-4] == 'G':
        return -4
    self_cmpl_result = check_max_self_complement_hybird(''.join(seq))
    if self_cmpl_result[0] >= Parameter.primer_self_max_cmpl_num or \
            self_cmpl_result[1] >= Parameter.primer_self_max_ctn_cmpl_num:
        return -5
    return 0

def create_random_prime(length, max_loop=100000, if_consider_Tm=True) -> str:
    """
    创建一个随机引物，使其满足基本引物约束条件
    :param length: 目标引物长度
    :return: 目标引物
    """
    satisfy = False
    loop_time = 0
    result = None

    a = 0
    b = 0
    c = 0
    d = 0
    e = 0

    while (not satisfy) and (loop_time < max_loop):
        loop_time += 1
        list_prime = list(create_random_strand(length-1, GC_content=0.50))
        last_base_index = random.randint(0, 2)
        if last_base_index == 0:
            list_prime.append('T')
        elif last_base_index == 1:
            list_prime.append('C')
        else:
            list_prime.append('G')
        #---------------------检查引物生化指标-----------------
        if lapi.get_GC_content(list_prime) < Parameter.primer_GC_content_min or\
                lapi.get_GC_content(list_prime) > Parameter.primer_GC_content_max:
            a+=1
            continue
        if if_consider_Tm:
            if abs(lapi.simple_get_Tm(''.join(list_prime)) - 57.5) >= Parameter.primer_Tm_bias_from_57_5_degree:
                b+=1
                continue
        if primer3.calcHairpin(''.join(list_prime)).structure_found is True:
            c+=1
            continue
        if list_prime[-2] == list_prime[-3] == list_prime[-4] == 'C' or\
            list_prime[-2] == list_prime[-3] == list_prime[-4] == 'G':
            d+=1
            continue
        self_cmpl_result = check_max_self_complement_hybird(''.join(list_prime))
        if self_cmpl_result[0] >= Parameter.primer_self_max_cmpl_num or\
                self_cmpl_result[1] >= Parameter.primer_self_max_ctn_cmpl_num:
            e+=1
            continue

        satisfy = True
        result = list_prime

    return ''.join(result)

def primer_library_statistic_use_only__create_random_primer(length, max_loop=100000, if_consider_Tm=True)->list:
    """
        仅针对primer library statistic文件使用的引物设计函数
        :param length: 目标引物长度
        :return: 目标引物
        """
    loop_time = 0
    result = 0
    pool = []

    while loop_time < max_loop:
        loop_time += 1
        list_prime = list(create_random_strand(length, GC_content=0.50))
        # ---------------------检查引物生化指标-----------------
        if list_prime[-1] == 'A':
            continue
        if lapi.get_GC_content(list_prime) < Parameter.primer_GC_content_min or \
                lapi.get_GC_content(list_prime) > Parameter.primer_GC_content_max:
            continue
        if if_consider_Tm:
            if abs(lapi.simple_get_Tm(''.join(list_prime)) - 57.5) >= Parameter.primer_Tm_bias_from_57_5_degree:
                continue
        if primer3.calcHairpin(''.join(list_prime)).structure_found is True:
            continue
        if list_prime[-2] == list_prime[-3] == list_prime[-4] == 'C' or \
                list_prime[-2] == list_prime[-3] == list_prime[-4] == 'G':
            continue
        self_cmpl_result = check_max_self_complement_hybird(''.join(list_prime))
        if self_cmpl_result[0] >= Parameter.primer_self_max_cmpl_num or \
                self_cmpl_result[1] >= Parameter.primer_self_max_ctn_cmpl_num:
            continue
        satisfy = True
        for p in pool:
            if p == ''.join(list_prime):
                satisfy = False
                break
            if lapi.check_hamming_distance_between_two_seqs(p, ''.join(list_prime)) <= Parameter.primer_lib_min_hamming_distance:
                satisfy = False
                break
        if satisfy:
            pool.append(''.join(list_prime))
        else:
            continue

        result += 1

    return pool

def check_max_self_complement_hybird(seq, show_structure=True) -> (int, int):
    """
    搜索单个序列自身的最大非特异性杂交情况
    :param seq: 序列
    :param show_structure: 是否显示杂交结构
    :return: 最大互补碱基数量, 最大连续互补碱基数
    """
    max_cmpl_num = 0
    max_ctn_cmpl_num = 0 # 最大连续互补碱基数
    output_structure = ''
    for i in range(len(seq)-1):
        l1 = list(seq)[:i+1]
        l1.reverse()
        l2 = list(seq)[i+1:]
        time = len(l1) if len(l1) < len(l2) else len(l2)
        cmpl_num = 0  # 互补碱基数
        cmpl_pos = []  # 互补碱基位置

        for j in range(time):
            if lapi.check_if_base_complement(l1[j], l2[j]):
                cmpl_num += 1
                cmpl_pos.append(j)
        # 最长连续数
        _t = lapi.check_max_continues_length(cmpl_pos)
        if max_ctn_cmpl_num < _t:
            max_ctn_cmpl_num = _t

        if cmpl_num > max_cmpl_num:
            max_cmpl_num = cmpl_num
            if show_structure:
                t1 = ''.join(l1)
                t2 = ''.join(l2)
                tm = []
                m = 0

                for k in range(time):
                    if m < len(cmpl_pos) and k == cmpl_pos[m]:
                        tm.append('|')
                        m += 1
                    else:
                        tm.append(' ')
                output_structure = t1 + '\n' + ''.join(tm) + '\n' + t2
    if show_structure:
        return max_cmpl_num, max_ctn_cmpl_num, output_structure
    return max_cmpl_num, max_ctn_cmpl_num

def check_max_2seqs_complement_hybird(seq1, seq2, show_structure=True) -> (int, int, float, int):
    """
    搜索两个序列的最大非特异性杂交情况,输入序列都从5‘ -> 3’
    :param seq1: 序列1
    :param seq2: 序列2
    :param show_structure: 是否显示杂交结构
    :return: 最大互补碱基数量，最大连续互补碱基数量, 引物与非特异性扩增区的序列最大同源性，引物3'末端最大连续互补碱基数
    """
    max_cmpl_num = 0
    max_ctn_cmpl_num = 0
    max_3_end_ctn_cmpl_num = 0

    time = abs(len(seq1) - len(seq2))
    l1 = list(seq1)
    l2 = list(seq2)
    l2.reverse()
    cmpl_num = 0  # 互补碱基数
    cmpl_pos = []  # 互补碱基位置
    output_structure = ''

    # if time == 0:
    #     for i in range(len(l1)):
    #         if lapi.check_if_base_complement(l1[i], l2[i]):
    #             max_cmpl_num += 1
    #             cmpl_pos.append(i)
    #
    #     max_ctn_cmpl_num = lapi.check_max_continues_length(cmpl_pos)
    #     if show_structure:
    #         t1 = ''.join(l1)
    #         t2 = ''.join(l2)
    #         tm = []
    #         m = 0
    #
    #         for k in range(len(l1)):
    #             if m < len(cmpl_pos) and k == cmpl_pos[m]:
    #                 tm.append('|')
    #                 m += 1
    #             else:
    #                 tm.append(' ')
    #         output_structure = t1 + '\n' + ''.join(tm) + '\n' + t2
    #         print(output_structure)
    #     print((max_cmpl_num, max_ctn_cmpl_num))
    #     return max_cmpl_num, max_ctn_cmpl_num
    #
    # else:
    short_seq = deepcopy(l1) if len(seq1) <= len(seq2) else deepcopy(l2)
    long_seq = deepcopy(l2) if len(short_seq) == len(seq1) else deepcopy(l1)

    max_cmpl_pos = []
    max_cmpl_gap = 0
    max_cmpl_long_seq_gap = 0

    max_ctn_cmpl_pos = []
    max_ctn_cmpl_gap = 0
    max_ctn_long_seq_gap = 0
    # 前端互补情况判定
    for i in range(len(short_seq)-1):

        cmpl_pos.clear()
        cmpl_num = 0
        ctn_cmpl_num = 0

        for j in range(i):
            if lapi.check_if_base_complement(short_seq[-i+j], long_seq[j]):
                cmpl_num += 1
                cmpl_pos.append(j)
            ctn_cmpl_num = lapi.check_max_continues_length(cmpl_pos)
            if max_cmpl_num < cmpl_num:
                max_cmpl_num = cmpl_num
                max_cmpl_pos = deepcopy(cmpl_pos)
                max_cmpl_long_seq_gap = len(short_seq) - i
            if max_ctn_cmpl_num < ctn_cmpl_num:
                max_ctn_cmpl_num = ctn_cmpl_num

            temp = deepcopy(cmpl_pos)
            for x in range(len(temp)):
                temp[x] += len(short_seq) - i
            end_3_ctn_cmpl_num = lapi.check_max_end_continues_length(temp, len(short_seq),
                                   Parameter.primer_3_end_max_ctn_cmpl_span)
            if max_3_end_ctn_cmpl_num < end_3_ctn_cmpl_num:
                max_3_end_ctn_cmpl_num = end_3_ctn_cmpl_num
                max_ctn_cmpl_pos = deepcopy(cmpl_pos)
                max_ctn_long_seq_gap = len(short_seq) - i


    # 中间互补情况判定
    for i in range(time):
        cmpl_pos.clear()
        cmpl_num = 0
        ctn_cmpl_num = 0

        for j in range(len(short_seq)):
            if lapi.check_if_base_complement(short_seq[j], long_seq[j+i]):
                cmpl_num += 1
                cmpl_pos.append(j+i)
        ctn_cmpl_num = lapi.check_max_continues_length(cmpl_pos)
        if max_cmpl_num < cmpl_num:
            max_cmpl_num = cmpl_num
            max_cmpl_pos = deepcopy(cmpl_pos)
            max_cmpl_gap = i
            max_cmpl_long_seq_gap = 0
        if max_ctn_cmpl_num < ctn_cmpl_num:
            max_ctn_cmpl_num = ctn_cmpl_num
            end_3_ctn_cmpl_num = lapi.check_max_end_continues_length(cmpl_pos, len(short_seq),
                              Parameter.primer_3_end_max_ctn_cmpl_span)

        temp = deepcopy(cmpl_pos)
        for x in range(len(temp)):
            temp[x] -= i
        end_3_ctn_cmpl_num = lapi.check_max_end_continues_length(temp, len(short_seq),
                                                                 Parameter.primer_3_end_max_ctn_cmpl_span)

        if max_3_end_ctn_cmpl_num < end_3_ctn_cmpl_num:
            max_3_end_ctn_cmpl_num = end_3_ctn_cmpl_num
            max_ctn_cmpl_pos = deepcopy(cmpl_pos)
            max_ctn_long_seq_gap = 0
            max_ctn_cmpl_gap = i

    # 后端互补情况判定
    for i in range(len(short_seq)):

        cmpl_pos.clear()
        cmpl_num = 0
        ctn_cmpl_num = 0

        for j in range(len(short_seq) - i):
            if len(long_seq)-len(short_seq)+i+j < len(long_seq) and\
            lapi.check_if_base_complement(short_seq[j], long_seq[len(long_seq)-len(short_seq)+i+j]):
                cmpl_num += 1
                cmpl_pos.append(len(long_seq)-len(short_seq)+i+j)
            ctn_cmpl_num = lapi.check_max_continues_length(cmpl_pos)
            if max_cmpl_num < cmpl_num:
                max_cmpl_num = cmpl_num
                max_cmpl_pos = deepcopy(cmpl_pos)
                max_cmpl_long_seq_gap = 0
                max_cmpl_gap = len(long_seq) - len(short_seq) + i
            if max_ctn_cmpl_num < ctn_cmpl_num:
                max_ctn_cmpl_num = ctn_cmpl_num
                end_3_ctn_cmpl_num = lapi.check_max_end_continues_length(cmpl_pos, len(short_seq),
                                            Parameter.primer_3_end_max_ctn_cmpl_span)

            temp = deepcopy(cmpl_pos)
            for x in range(len(temp)):
                temp[x] -= len(long_seq) - len(short_seq) + i
            end_3_ctn_cmpl_num = lapi.check_max_end_continues_length(temp, len(short_seq),
                                                                     Parameter.primer_3_end_max_ctn_cmpl_span)

            if max_3_end_ctn_cmpl_num < end_3_ctn_cmpl_num:
                max_3_end_ctn_cmpl_num = end_3_ctn_cmpl_num
                max_ctn_long_seq_gap = 0
                max_ctn_cmpl_pos = deepcopy(cmpl_pos)
                max_ctn_cmpl_gap = len(long_seq) - len(short_seq) + i



    if show_structure:
        l1 = deepcopy(list(short_seq))
        l2 = deepcopy(list(long_seq))
        for _ in range(max_cmpl_gap):
            l1.insert(0, ' ')
        for _ in range(max_cmpl_long_seq_gap):
            l2.insert(0, ' ')
        t1 = ''.join(l1)
        t2 = ''.join(l2)
        tm = [' ' for _ in range(max_cmpl_gap + max_cmpl_long_seq_gap)]
        m = 0

        for k in range(len(short_seq) - max_cmpl_long_seq_gap):
            if m < len(max_cmpl_pos) and k + max_cmpl_gap == max_cmpl_pos[m]:
                tm.append('|')
                m += 1
            else:
                tm.append(' ')
        output_structure = t1 + '\n' + ''.join(tm) + '\n' + t2
        print(output_structure)


    # print((max_cmpl_num, max_ctn_cmpl_num, max_cmpl_num / len(short_seq), max_3_end_ctn_cmpl_num))
    return max_cmpl_num, max_ctn_cmpl_num, max_cmpl_num / len(short_seq), max_3_end_ctn_cmpl_num

def design_primer_library(size, primer_length, max_loop=100000, output_primer_number_realtime=True, output_process=True, consider_tm=True) -> list:
    """
    设计引物库
    :param size: 引物库大小
    :param primer_length: 引物长度
    :param output_primer_number_realtime: 是否实时输出已经设计出的引物数量
    :param output_process: 是否输出遍历百分比
    :return:
    """
    primer_list = []

    primer_list.append(create_random_prime(primer_length, 100000, consider_tm))
    loop_time = 0
    while len(primer_list) < size and loop_time < max_loop:
        loop_time += 1
        if primer_length < 18:
            t_primer = create_random_prime(primer_length, 100000, if_consider_Tm=False)
        else:
            t_primer = create_random_prime(primer_length, 100000)
        satisfy = True
        for primer in primer_list:
            check_result = check_max_2seqs_complement_hybird(primer, t_primer, show_structure=False)
            if check_result[0] >= Parameter.primer_lib_max_cmpl_ratio*\
            primer_length or check_result[1] >= Parameter.primer_lib_max_ctn_cmpl_num:
                satisfy = False
                break
            if lapi.check_hamming_distance_between_two_seqs(primer, t_primer) <= Parameter.primer_lib_min_hamming_distance:
                satisfy = False
                break
        if satisfy:
            primer_list.append(t_primer)
            # loop_time = 0
            if output_primer_number_realtime:
                print("成功设计出{}个{}-mer引物".format(len(primer_list), len(primer)))
        if int(loop_time % (max_loop / 10)) == 0:
            if output_process:
                print("已完成-{}%".format(loop_time*100 / max_loop), end=' ')

    return primer_list

def pick_final_primer_library(meta_primer_lib) ->list:
    """
    挑选最终满足正交性约束的引物
    :param meta_primer_lib: 原始引物库
    :return: 满足正交性约束的引物库
    """
    pb = ProgressBar()
    final_primer_lib = []
    try:
        meta_primer_lib = list(meta_primer_lib)
        primer_length = len(meta_primer_lib[0])

        final_primer_lib.append(meta_primer_lib.pop(0))
    except IndexError:
        meta_primer_lib = []
    pb.currval = 0
    pb.maxval = len(meta_primer_lib)
    pb.start()
    for p in meta_primer_lib:
        satisfy = True
        for fp in final_primer_lib:
            check_result = check_max_2seqs_complement_hybird(p, fp, False)
            if check_result[0] >= Parameter.primer_lib_max_cmpl_ratio * \
                    primer_length or check_result[1] >= Parameter.primer_lib_max_ctn_cmpl_num:
                satisfy = False
                break
        if satisfy:
            final_primer_lib.append(p)

        pb.update(meta_primer_lib.index(p)+1)
    pb.finish()
    return final_primer_lib

def assemble_primer_and_payload(forward_primer_list, payload, reversed_primer_list) -> (str, str):
    """
    组装正向引物、有效载荷和反向引物
    :param forward_primer_list: 正向引物列表
    :param payload: 有效载荷
    :param reversed_primer_list:反向引物列表
    :return: 完整序列
    """
    result = ""
    pointer = []
    for primer in forward_primer_list:
        result += primer
        for p in range(len(primer)):
            if p != len(primer) - 1:
                pointer.append(' ')
            else:
                pointer.append('^')
    result += payload
    for p in payload:
        pointer.append('_')

    for primer in reversed_primer_list:
        result += get_complement_strand(primer)
        for p in range(len(primer)):
            if p != len(primer) - 1:
                pointer.append(' ')
            else:
                pointer.append('^')
    return result, ''.join(pointer)
