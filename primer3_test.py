import random

import primer3
import pandas
import Functions.Strand as strand
import Functions.Primer3Tools as primer3tools
import threading
import time
from Functions.SequenceIterator import SeqGenerator
import Functions.Encoding as Encoding

#strand.create_random_prime(20, 5000)




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
#
# universal_forward_primer = 'AGTGAGGGTTCGTGGGTGGT'
# universal_reversed_primer = 'GCTTGAGTTGAGTGGGGAGT'
# forward_primers = ['TCGGGTCTACGGCCAGATTT', 'GGAGAGTGAGGAGGGCTTCT', 'TGAGATTTCAGGGGGCGGCT']
# reversed_primers = ['TCGGGCTTGTGGTGTGGGAT', 'GAGAGGTGGAGTTGCGAGTT', 'TGTGGGTGAGGTGTGAGAGT']
# features = ['GAGGGGTGGGGGGTTTATTTGGGGGTGGATGGTGGTTATT',
#             'GCGTGTGTTGGAGGTCGTTTACGGATTGGGGCTGCGATGT',
#             'GGCGTGAGGGAGGTTGTTGTTTGTGCGGTCGGGGGTGTTT',
#             'TTGGGGTGGTTTTGGTGCCTTGAGTGTGCGGTGTGGATGT',
#             'CTTTTGGTGGGGGGTCGTTTGCGGAGTTAGGAGGGGTTTT',
#             'GGGCGTTGTTGTGCTGCGTTGGTGCTTATTAGGGGGGGTT',
#             'TGTGTTCTTTCGGCGGGTGTCGGTGAGTGGTTGTGGGGAT',
#             'CTTGCGGGGGGAGAGGTTTTTGGGTCGGGTGAGGTTTGTT',
#             'AGAGTGTGGGTTTTGGGGCTTGGTGGTCTTGGGAGGTTGT',
#             'TAGGAGTGGGGTTTGCGGCTGTGTTTTTGTGCGGGGGAGT',
#             'GGTGGGGAGAGAGTTCGTTTTGCGGGGATGGTGGGATTCT',
#             'TGTGGGTGGTTGAGGGGTATATTGTGGGGGGGTTCGTTCT']
#
# l1 = []
# l1.append(universal_forward_primer)
# l1.append(forward_primers[0])
# l1.append(forward_primers[1])
# l2 = []
# l2.append(reversed_primers[2])
# l2.append(reversed_primers[1])
# l2.append(universal_reversed_primer)
#
# s = strand.assemble_primer_and_payload(l1, features[9], l2)
# print(s[0])
# print(s[1])




# with open('20merPrimer.txt', 'r') as f:
#     content = f.read()
#     content = content.replace('[', '')
#     content = content.replace(']', '')
#     content = content.replace(' ', '')
#     content = content.replace('\'', '')
#     content = content.split(',')
# for item in content:
#     item = strand.get_complement_strand(item, False)
#     item = item[:-1] + 'T'
#     print(strand.check_if_satisfy_primer_principle(item))


# seqGenerator = SeqGenerator(20, 2**38)
# seqIter = iter(seqGenerator)
# print(next(seqIter))
# times = time.time()
# for i in range(2**40):
#     seq = next(seqIter)
#     strand.check_if_satisfy_primer_principle(seq)
#     if i % (2**20) == 0:
#         print(time.time()-times)
# print(time.time()-times)

f = open('FinalPrimer.txt', 'r')
content = f.read()
content = content.replace('[', '')
content = content.replace(']', '')
content = content.replace('\'', '')
content = content.replace(' ', '')

content = content.split(',')

shield_bases = 'ATC'
UFP = content[0]
URP = content[1]
FPs = content[2:5]
RPs = content[5:8]

forward_primer_list = []
reversed_primer_list = []

for fp1 in FPs:
    for fp2 in FPs:
        forward_primer_list.append(fp1 + fp2)
for rp1 in RPs:
    for rp2 in RPs:
        reversed_primer_list.append(rp1 + rp2)

print("通用上游引物:", UFP)
print("通用下游引物:", URP)
print("正向引物库:", FPs)
print("反向引物库:", RPs)
payloads = []
# for _ in range(81):
#     payloads.append(strand.create_random_strand(60))
word_list = [
    "Hello World!!!!", "I Love TJU(^ ^)", "uncopyrightable", "abiogenetically", "abnormalization",
    "easternizations", "abiogenetically", "aboriginalities", "absorbabilities", "spectroscopists",
    "overrefinements", "inexplicability", "heartsicknesses", "internationally", "oligonucleotide",
    "poststimulatory", "jurisprudential", "valetudinarians", "reformabilities", "interiorization",
    "overencouraging", "platitudinously", "transsexualisms", "nonimplications", "lightsomenesses",
    "infundibuliform", "deliriousnesses", "geotectonically", "hyperaggressive", "palatalizations",
    "unimaginatively", "alkalinizations", "denationalizing", "nonfulfillments", "fraternizations",
    "hypnotherapists", "intraperitoneal", "intrapreneurial", "sorrowfulnesses", "noncompressible",
    "hyperparasitism", "interventionist", "postapocalyptic", "industriousness", "nonexploitation",
    "nonbiographical", "cinematographer", "protohistorians", "photoelectronic", "gentlemanliness",
    "morphologically", "microporosities", "sensationalists", "intramuscularly", "monopolizations",
    "phototoxicities", "equalitarianism", "unobjectionable", "unreliabilities", "seaworthinesses",
    "ariboflavinoses", "infusiblenesses", "disappointingly", "disinterestedly", "psychiatrically",
    "undemonstrative", "incorruptnesses", "proletarianised", "preimplantation", "circularization",
    "superstimulated", "micropublishers", "straightforward", "overaccentuated", "unselfishnesses",
    "subdevelopments", "renationalizing", "lignocellulosic", "lexicalisations", "quatercentenary"
]

for word in word_list:
    payloads.append(Encoding.encode_string_to_dna_seq(word))

final_strands = []
t = 0
for i in range(len(forward_primer_list)):
    for j in range(len(reversed_primer_list)):
        if forward_primer_list[i][:20] == forward_primer_list[i][20:] and \
            reversed_primer_list[j][:20] == reversed_primer_list[j][20:]:
            continue
        payload = payloads[t]
        t += 1
        std = strand.assemble_primer_and_payload([UFP, forward_primer_list[i]], payload, [reversed_primer_list[j], URP])[0]
        std = shield_bases + std + shield_bases
        final_strands.append(str(i)+std+str(j))

for s in final_strands:
    print(s, len(s)-2, "File-"+str(final_strands.index(s)))

print(len(final_strands))

"""
正向引物库: ['GAAGGTCTGTCATGGTTCTG', 'TGGGACTGATGTGGACTGCC', 'GGTCATAGCCTCATCCAGCG']
反向引物库: ['TTGGCACTGAATCACGCGTC', 'AGTTCTAACGGCAGGATAGT', 'TTTTTAATCCCCTCCATCCG']

带接头正向引物库: ['TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGAAGGTCTGTCATGGTTCTG',
            'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTGGGACTGATGTGGACTGCC',
            'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGGTCATAGCCTCATCCAGCG']
带接头反向引物库: ['GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTTGGCACTGAATCACGCGTC',
            'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGAGTTCTAACGGCAGGATAGT',
            'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTTTTTAATCCCCTCCATCCG']
"""