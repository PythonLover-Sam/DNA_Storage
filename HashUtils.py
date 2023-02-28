def bkdrHash(input: str):
    seed = 131
    hash = 0
    for char in input:
        hash = hash * seed + ord(char)
        hash &= 0xFFFFFFFF
    return hash & 0x7FFFFFFF
