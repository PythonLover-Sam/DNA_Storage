import pandas


"""
将Primer3 引物设计结果字典转换为pandas中的Dataframe格式
"""
def primer_design_result_to_pretty_dataFrame(result):
    result_table = {} 
    keys = []

    for key in result:
        if '0' in key:
            keys.append(key.replace('_0', ''))

    for i in range(result['PRIMER_PAIR_NUM_RETURNED']):
        new_key = "Primer_" + str(i)
        new_list = []
        for key in result:
            if str(i) in key:
                new_list.append(result[key])
        result_table[new_key] = new_list

    df = pandas.DataFrame(result_table, index=keys)
    return df

