import random

import primer3
import pandas
import Functions.Strand as strand
import Functions.Primer3Tools as primer3tools
import threading
import time
import Functions.LowLevelAPI as lapi
from Functions.SequenceIterator import SeqGenerator
import Parameters


output_file = open('output.txt', 'a')

def search_every_task(length, info):
    print(info)
    pass_list = []
    seqGenerator = SeqGenerator(length)
    seqIter = iter(seqGenerator)

    _pass = 0
    _fail = 0

    for i in range(3*4**(length-1)):
        seq = next(seqIter)
        result = strand.check_if_satisfy_primer_principle(seq, False)
        if result == 0:
            satisfy = True
            for p in pass_list:
                if lapi.check_hamming_distance_between_two_seqs(p, seq) <= Parameters.Parameter.primer_lib_min_hamming_distance:
                    satisfy = False
                    break
            if satisfy:
                pass_list.append(seq)
                _pass += 1
            else:
                _fail += 1


        else:
            _fail += 1
    print(info+' 完成遍历 引物长度[{}] 满足的引物数[{}] 不满足的引物数[{}] 通过率[{}%]'.format(length, _pass, _fail+4**(length-1), _pass*100/(4**length)))
    output_file.write(info + ' 引物长度[{}] 满足的引物数[{}] 不满足的引物数[{}] 通过率[{}%] \n'.format(length, _pass, _fail+4**(length-1), _pass*100/(4**length)))

def search_random_task(length, pool_size, info):
    print(info+ 'start')
    #result = strand.design_primer_library(9999999999999, length, pool_size, False, True, False)
    result = len(strand.primer_library_statistic_use_only__create_random_primer(length, pool_size, False))
    print(info + ' 完成随机引物设计 引物长度[{}] 满足的引物数[{}] 不满足的引物数[{}] 通过率[{}%]'.format(
        length, result, pool_size-result, result*100/pool_size))
    output_file.write(info + ' ,{} ,{}, {}, {}, {}% \n'.format(
        length, pool_size, result, pool_size-result, result*100/pool_size))

ts = []

# ts.append(threading.Thread(target=search_every_task, args=(5, 'ID1')))
# ts.append(threading.Thread(target=search_every_task, args=(6, 'ID2')))
# ts.append(threading.Thread(target=search_every_task, args=(7, 'ID3')))
# ts.append(threading.Thread(target=search_every_task, args=(8, 'ID4')))
# ts.append(threading.Thread(target=search_every_task, args=(9, 'ID5')))
# ts.append(threading.Thread(target=search_every_task, args=(10, 'ID6')))
# ts.append(threading.Thread(target=search_every_task, args=(11, 'ID7')))
# ts.append(threading.Thread(target=search_every_task, args=(12, 'ID8')))



# pool_size = 1e5
# for i in range(5):
#     p_length = 26 + i
#     info = str(pool_size)+', ' +str(p_length)
#     for j in range(10):
#         ts.append(threading.Thread(target=search_random_task, args=(p_length, pool_size, info)))
#
#
#
# for t in ts:
#     t.start()

# def design_20mer_primer():
#
#     print(strand.design_primer_library(1000, 20, max_loop=300000))
# def design_25mer_primer():
#
#     print(strand.design_primer_library(1000, 25, max_loop=300000))
#
# t1 = threading.Thread(target=design_20mer_primer)
# t2 = threading.Thread(target=design_25mer_primer)
#
# t1.start()
# t2.start()

number = 0
for i in range(10):

    result = strand.primer_library_statistic_use_only__create_random_primer(25, 1e4)
    result = strand.pick_final_primer_library(result)
    number += len(result)
number /= 10
f = open('20mer__primerLibrary.txt', 'a')
f.write('1e3, {}\n'.format(number))
f.close()