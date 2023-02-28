class Parameter:

    primer_GC_content_min = 0.45            # 引物最小GC含量
    primer_GC_content_max = 0.60            # 引物最大GC含量
    primer_Tm_bias_from_57_5_degree = 25     # 引物Tm值偏离57.5摄氏度的最大程度
    primer_self_max_cmpl_num = 5            # 引物自身最大互补碱基数(不包含)
    primer_self_max_ctn_cmpl_num = 4        # 引物自身最大连续互补碱基数(不包含)
    primer_3_end_max_ctn_cmpl_num = 8       # 引物3‘端最大互补碱基数
    primer_3_end_max_ctn_cmpl_span = 2      # 引物3’端最大互补碱基数考虑跨度
    primer_lib_max_cmpl_ratio = 0.3         # 引物之间最大同源性
    primer_lib_max_ctn_cmpl_num = 5         # 引物之间最大连续互补碱基数
    primer_lib_min_hamming_distance = 5     # 引物之间最小汉明距离