def decode_dna_seq_to_string(dna_seq, force_decode=False):
    """

    :param dna_seq:
    :param force_decode: 是否强制解码，不考虑奇偶校验
    :return:
    """
    result = ""
    bin_seq = ""
    for base in dna_seq:
        bin_seq += get_dna_base_to_bin(base)
    if not force_decode:
        check_code = bin_seq[-2:]
        num0 = 0
        num1 = 0
        for bin_code in bin_seq[:-2]:
            if bin_code == '0':
                num0 += 1
            else:
                num1 += 1
        if (num0 + (int(check_code[0])+1)%2)%2 != 1 or (num1 + int(check_code[1]))%2 != 1:
            print((num0 + (int(check_code[0])+1)%2)%2, (num1 + int(check_code[1]))%2)
            raise Exception("校验出错")
    for i in range(len(bin_seq) // 8):
        char = chr(int(bin_seq[i*8:i*8+8], 2))
        result += char
    return result
def encode_string_to_dna_seq(string):
    result = ""
    bin_seq = ""
    bin_seq += get_string_ascii_bin(string)
    bin_seq += get_check_code(bin_seq)
    if len(bin_seq) % 2 != 0:
        raise ValueError("二进制序列长度不是偶数")
    for i in range(len(bin_seq)//2):
        result += get_bin_to_dna_base(bin_seq[2*i: 2*i+2])
    return result
def get_dna_base_to_bin(dna_base):
    if dna_base == 'A':
        return '00'
    elif dna_base == 'C':
        return '01'
    elif dna_base == 'T':
        return '10'
    else:
        return '11'
def get_bin_to_dna_base(str_bin):
    if str_bin == '00':
        return 'A'
    elif str_bin == '01':
        return 'C'
    elif str_bin == '10':
        return 'T'
    else:
        return 'G'
def get_check_code(ascii_seq):  # 奇校验
    num0 = 0
    num1 = 0
    for char in ascii_seq:
        if char == '0':
            num0+=1
        else:
            num1+=1
    result = ''
    if num0 % 2 == 0:
        result += '0'
    else:
        result += '1'
    if num1 % 2 == 0:
        result += '1'
    else:
        result += '0'
    return result
def get_string_ascii_bin(string):
    result = ""
    for char in string:
        result += get_char_ascii_bin(char)
    return result
def get_char_ascii_bin(char):
    result = []
    charNum = ord(char)
    for _ in range(8):
        result.append(str(int(charNum % 2)))
        charNum /= 2
    result.reverse()
    return ''.join(result)


if __name__ == "__main__":
    print(1)
