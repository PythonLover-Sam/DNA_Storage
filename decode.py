from Functions import Encoding as encoding


seq = "CATACTCCCTGACTGACTGGATAACCCGCTGGCGATCTGACTCAATACATACATACATACT"

print(encoding.decode_dna_seq_to_string(seq))

seq = "CATCATAACAGACTGGCGCTCTCCATAACCCACATTCCCCATTACCGTATAACCGTATTCT"

print(encoding.decode_dna_seq_to_string(seq))
seq = "CTACCTATCGAGCTGGCGATCTATCTACCTATCTTCCTGACTTCCGCACTTCCTCCCGAGT"
print(encoding.decode_dna_seq_to_string(seq))

"""
UFP--FP1--FP1--I Love TJU(^ ^)--RP3--RP1--URP


UFP--FP1--FP2--absorbabilities--RP3--RP1--URP

"""
