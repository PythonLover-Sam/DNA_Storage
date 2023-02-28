# 序列迭代器 遍历全部序列
def base_index_to_base(base_index):
    if base_index == 0:
        return 'A'
    elif base_index == 1:
        return 'C'
    elif base_index == 2:
        return 'G'
    elif base_index == 3:
        return 'T'
    else:
        raise ValueError("碱基数字表示不正确，需要范围0-3，但是得到的是{}".format(base_index))


class SeqGenerator:
    def __init__(self, length, seq_index=-1):
        if seq_index == -1:
            self.seq_index = 2**(2*length-2)
        else:
            self.seq_index = seq_index
        self.length = length
        self.seq = self.translate()
    def __iter__(self):
        return self
    def __next__(self):
        seq = self.seq
        self.seq_index += 1
        self.seq = self.translate()
        return seq

    def translate(self):
        seq = []
        _seq_index = self.seq_index
        for i in range(self.length):
            base_index = (_seq_index | 0b00) & 0b11
            seq.append(base_index_to_base(base_index))
            _seq_index >>= 2
        return ''.join(seq)
