import Functions.Strand as strand
import Functions.LowLevelAPI as lapi
from Functions.SequenceIterator import SeqGenerator
import Parameters
import getopt, sys
import progressbar

p = progressbar.ProgressBar()
loop_time = 100000
single_loop_time = 1e4
base_num = 20
opts,args = getopt.getopt(sys.argv[1:],'-h-f:-v-t:-s:-n:',['help','filename=','version', 'loop_time=', 'single_loop_time=', 'base_num='])
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print("[*] Help info")
        exit()
    if opt_name in ('-v','--version'):
        print("[*] Version is 0.01 ")
        exit()
    if opt_name in ('-f','--filename'):
        fileName = opt_value
        print("[*] Filename is ",fileName)
        # do something
    if opt_name in('-t', '--loop_time'):
        loop_time = int(opt_value)
    if opt_name in ('-s', '--single_loop_time'):
        single_loop_time = int(opt_value)
    if opt_name in ('-n', '--base_num'):
        base_num = int(opt_value)

number = 0
for i in range(loop_time):
    print(i)
    result = strand.primer_library_statistic_use_only__create_random_primer(base_num, single_loop_time)
    result = strand.pick_final_primer_library(result)
    if loop_time == 1:
        f = open(fileName, 'a')
        f.write(str(result))
        f.close()
    number += len(result)
number /= loop_time
f = open('20mer__primerLibrary.txt', 'a')
f.write('{}, {}, {}\n'.format(base_num, single_loop_time, number))
f.close()