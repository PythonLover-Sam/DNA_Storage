# import Functions.Encoding as encoding
# import Functions.Strand as strand
# import Functions.LowLevelAPI as lapi
#
#
# FPs = ["AGAATCGCTGTGCTTCTTTG",  # 原始正向引物
#        "AGAATCGCTGTGCTTCTCTG",  # 1碱基差异  right
#        #                 ^
#        "AGAATTGCTGTGCATCTCTG",  # 2碱基差异 middle
#        #     ^       ^
#        "ATAATCGCTGTGCTTCTGTG",  # 2碱基差异 side
#        # ^               ^
#        "AGAATTGCTCTGCATCTCTG",  # 3碱基差异 middle
#        #     ^   ^   ^
#        "ACAATCGATGTGCTACTTAG",  # 4碱基差异 uniform
#        # ^     ^      ^   ^
#        "AGAAGAGTTGCACTTCTTTG",  # 5碱基差异 middle
#        #    ^^ ^  ^^
#        "AGAATCGCTGTGCTTCTTTG",  # 原始正向引物
#        "ATAAGCGCTATGGTTGTTCG",  # 6碱基差异
#        # ^  ^    ^  ^  ^  ^
#        "AGTATCTATATGGCTCATCG",  # 8碱基差异
#        #  ^   ^^ ^  ^^  ^ ^
#        ]
# RPs = ["GATTGCAGCAAGAACTCCGT",  # 原始反向引物
#        "GAGTGCAGCAAGAACTCCGT",  # 1碱基差异  left
#        #  ^
#        "AAGTGCAGCAAGAACTCCGC",  # 2碱基差异  side
#        #^                  ^
#        "GATTGCAGCAAGAACTCCGT",  # 原始反向引物 无差异
#        "TATTGCAGCCAGAACTCCGG",  # 3碱基差异  uniform
#        #^        ^         ^
#        "GATTGTAACATGAATTCCGT",  # 4碱基差异  middle
#        #     ^ ^  ^   ^
#        "GATTGCAGCAAGAACTCCGT",  # 原始反向引物 无差异
#        "GGTTGCGGCAAAAACTTCGC",  # 5碱基差异  uniform
#        # ^    ^    ^    ^  ^
#        "GAGTGCTACAATAAAGCCGT",  # 6碱基差异
#        #  ^   ^^   ^  ^^
#        "GATTGCAGCAAGAACTCCGT",  # 原始反向引物 无差异
#        ]
#
# shield = "ATC"
#
# sentences = [
#     "this is the first original strand 0001",
#     "ok very good performance Geforce: 0002",
#     "3x faster than previous laptops : 0003",
#     "Nvidia GeForce 3090Ti good GPU *: 0004",
#     "Xiao Chun is the apple of my eye: 0005",
#     "this is the other one strand num: 0006",
#     "this is the other one strand num: 0007",
#     "this is the other one strand num: 0008",
#     "this is the other one strand num: 0009",
#     "this is the other one strand num: 0010",
# ]
#
# payloads = []
#
# for sentence in sentences:
#     print(len(encoding.encode_string_to_dna_seq(sentence)))
#     payloads.append(encoding.encode_string_to_dna_seq(sentence))
#
# print(payloads)
#
# for payload in payloads:
#     print(encoding.decode_dna_seq_to_string(payload))
#
# final_strands = []
# for fp, payload, rp in zip(FPs, payloads, RPs):
#     result = strand.assemble_primer_and_payload([fp], payload, [rp])[0]
#     result = shield + result + shield
#     print(result, len(result))
#

import Functions.Strand as strand
import Functions.Encoding as encoding
# print(strand.get_complement_strand("CTTTCGCCCGATCTTCCGAGCGAACGATCGCCCTCACTCCCTGTCGCACTTCCTACCTGAT"))

seq = encoding.encode_string_to_dna_seq("")
print(encoding.decode_dna_seq_to_string(seq[:], True))

