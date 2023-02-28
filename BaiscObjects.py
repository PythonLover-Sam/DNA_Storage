from HashPrimerIndexMatrix import HashUtils


class HashPrimerItem:

    pattern = None

    def __init__(self):
        if HashPrimerItem.pattern is None:
            raise ValueError("必须在创建哈希引物元素之前设置HashPrimerItem.pattern类变量的初始化值，需要传入一个元组描述矩阵形状")
        self.primer = ""
        self.itemIndex = 0
        self.file = None
        self.coordinate = ()

    def getPrimer(self):
        return self.primer

    def getPrimerLength(self):
        return len(self.primer)

    def setPrimer(self, primer: str):
        self.primer = primer
        return self

    def getItemIndex(self):
        return self.itemIndex

    def setItemIndex(self, itemIndex: int):
        self.itemIndex = itemIndex
        return self

    def getFile(self):
        return self.file

    def setFile(self, file: object):
        self.file = file
        return file

    def getCoordinate(self):
        return self.coordinate

    def setCoordinate(self, coor: tuple):
        self.coordinate = coor
        return self

    def getGroundTruthCoordinate(self) -> tuple:
        """
        获取通过引物index换算得到的本原坐标
        :return:
        """
        return (self.itemIndex % HashPrimerItem.pattern[0],
                self.itemIndex // HashPrimerItem.pattern[0] % (HashPrimerItem.pattern[1]),
                self.itemIndex // (HashPrimerItem.pattern[0] * HashPrimerItem.pattern[1]))

    def getCoordinateToColorMapIndex(self)-> int:
        """
        通过元素的实际坐标获取其伪itemIndex用于绘制3D图展示确定颜色映射表
        :return:
        """
        x, y, z = HashPrimerItem.pattern[0], HashPrimerItem.pattern[1], HashPrimerItem.pattern[2]
        return self.coordinate[0] * y*z + z * self.coordinate[1] + self.coordinate[2]

    def checkIfInOriginalCoordinate(self):
        """
        检查该元素是否在本来的位置
        :return:
        """
        return self.coordinate == self.getGroundTruthCoordinate()

class HashPrimerMatrix:

    def __init__(self, pattern: tuple, use_infinite_virtual_primer: bool):
        """
        哈希引物矩阵的构造函数
        :param pattern: 引物矩阵的形状
        :param use_infinite_virtual_primer: 是否使用无限虚拟引物
        """
        self.items: list[HashPrimerItem] = []
        self.pattern = pattern
        HashPrimerItem.pattern = pattern

        with open("20merPrimerLibraryFinal.txt", 'r') as f:
            p: list = eval(f.read())
            p.sort(key=lambda primer: PrimerToIndex(primer))

        for i in range(pattern[0] * pattern[1] * pattern[2]):
            item = HashPrimerItem().setItemIndex(i)
            item.setCoordinate(item.getGroundTruthCoordinate())
            if use_infinite_virtual_primer:
                item.setPrimer(p[i%len(p)])
            else:
                item.setPrimer(p[i])
            self.items.append(item)

    def coordinateToIndex(self, coordinate: tuple):
        """
        将给定的坐标元组转换为对应的此位置的本原索引值
        :param coordinate:
        :return:
        """
        x, y, z = self.pattern[0], self.pattern[1], self.pattern[2]
        return coordinate[0] + coordinate[1] * x + coordinate[2] * (x*y)

    def indexToCoordinate(self, index: int)-> tuple:
        """
        将给定的索引值转换为对应的此位置的本原坐标元组
        :param index:
        :return:
        """
        x, y, z = self.pattern[0], self.pattern[1], self.pattern[2]
        return index % x, index // x % y, index // (x * y)

    def trySaveFile(self, fileName: str, res: int, series: int):
        """
        尝试保存文件到item中
        :param fileName: 文件名（可哈希）
        :param res: 分辨率 整数序号代表
        :param series: 系列 整数序号代表
        :return:
        """
        hash = HashUtils.bkdrHash(fileName)
        x_pos = hash % self.pattern[0]
        original_pos = (x_pos, series, res)
        pos_item = self.items[self.coordinateToIndex(original_pos)]

        if pos_item.file is None:
            pos_item.file = fileName
            return self
        else:
            # 如果此处发生哈希冲突, 先延申探测，后线性探测
            check_pos = pos_item.coordinate
            last_item = pos_item
            i = 1
            while check_pos != original_pos:
                # 当没有遍历到碰撞末尾的时候
                last_item = self.items[self.coordinateToIndex(check_pos)]
                check_pos = last_item.coordinate
                i += 1
            # 如果遍历到末尾,新建一块存储空间，并且交换与前一个
            dy = 1
            while last_item.getGroundTruthCoordinate()[1] + dy < self.pattern[1]:

                # 如果沿着y轴还有空间
                new_pos_item = self.items[self.coordinateToIndex((
                    last_item.getGroundTruthCoordinate()[0],
                    last_item.getGroundTruthCoordinate()[1] + dy,
                    last_item.getGroundTruthCoordinate()[2]
                ))]
                if new_pos_item.file is None:
                    new_pos_item.file = fileName
                    temp = new_pos_item.coordinate
                    new_pos_item.coordinate = last_item.coordinate
                    last_item.coordinate = temp
                    return self
                else:
                    dy += 1
            # 如果沿着y轴已经没有空间，利用x轴向的线性探测法
            for i in range(self.pattern[0]):
                pos = self.items[self.coordinateToIndex(((x_pos + i) % self.pattern[0],
                                                        0,
                                                        last_item.getGroundTruthCoordinate()[2]))]
                if pos.file is None:
                    pos.file = fileName
                    temp = pos.coordinate
                    pos.coordinate = last_item.coordinate
                    last_item.coordinate = temp
                    return self
        return self
        # raise OverflowError("找不到可用存储空间")

def PrimerToIndex(primer: str)-> int:
    bins = ""
    for base in primer:
        bins += BaseToBin(base)
    return int(bins, 2)

def BaseToBin(base: str):
    if base == "A": return "00"
    if base == "C": return "01"
    if base == "G": return "10"
    if base == "T": return "11"

if __name__ == "__main__":
    HashPrimerItem.pattern = (2, 3, 3)
    a = HashPrimerItem().setItemIndex(21).getGroundTruthCoordinate()
    print(a)

    hashPrimerMatrix = HashPrimerMatrix((2, 3, 3))
    print(hashPrimerMatrix.items[3].checkIfInOriginalCoordinate())



