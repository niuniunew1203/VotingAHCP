# ==========================================
#  MACHINE LEARNING ANALYSIS
# ==========================================
import pandas as pd
import numpy as np
from Bio import SeqIO
import os
import sys
from itertools import product
import re
import math
def get_korder_base(chars, k):
    korder_chars = []
    chars_len = len(chars)
    for i in range(0, chars_len ** k):
        n = i
        bases = ''
        for j in range(0, k):
            base = chars[n % chars_len]
            n = n // chars_len
            bases += base
        korder_chars.append(bases[::-1])
    return korder_chars


def build_prot_mat(seq, str, k):
    seq_chars = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    str_chars = ['C', 'H', 'E']
    seq = re.sub(r"[UZOB]", "X", seq)
    seq_korder = get_korder_base(seq_chars, k)  # 20*20=400
    str_korder = get_korder_base(str_chars, k)
    row = len(seq_chars) ** k
    col = len(str_chars) ** k
    prot_fea_mat = np.zeros(shape=(row, col))

    for i in range(0, len(seq) - k + 1):
        s, e = i, i + k
        seq_base = seq[s:e]
        str_base = str[s:e]
        if 'X' in seq_base:
            prot_fea_mat[:][str_korder.index(str_base)] += 1 / len(seq_korder)
        else:
            row_index = seq_korder.index(seq_base)
            col_index = str_korder.index(str_base)
            prot_fea_mat[row_index][col_index] += 1
    return prot_fea_mat


def load_sequences(file_path):
    seq_dict = {}
    with open(file_path, 'r') as rf:
        seq = ''
        for line in rf:
            line = line.strip()
            if line[0] == '>':
                name = line[1:]
            else:
                seq = line.upper()
                seq_dict[name] = seq
    return seq_dict


letters = list('ACDEFGHIKLMNPQRSTVWY')


def readAAT(file):  # read AAT features from the AAT textfile
    try:
        aatdic = {}
        aatdata = open(file, 'r')
        for l in aatdata.readlines():
            aatdic[l.split()[0][0:3]] = float(l.split()[1])
        aatdata.close()
        return aatdic
    except:
        print("Error in reading AAT feature file. Please make sure that the AAT file is correctly formatted")
        sys.exit()


def readAAP(file):  # read AAP features from the AAP textfile
    try:
        aapdic = {}
        aapdata = open(file, 'r')
        for l in aapdata.readlines():
            aapdic[l.split()[0]] = float(l.split()[1])
        aapdata.close()
        return aapdic
    except:
        print("Error in reading AAP feature file. Please make sure that the AAP file is correctly formatted")
        sys.exit()


Amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
               'R', 'S', 'T', 'V', 'W', 'Y']
Amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
               'R', 'S', 'T', 'V', 'W', 'Y']
Amino_acids_ = list(product(Amino_acids, Amino_acids))
Amino_acids_ = [i[0] + i[1] for i in Amino_acids_]

# Amino_acids2 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
#                'R', 'S', 'T', 'V', 'W', 'Y','X']
# Amino_acids2_ = list(product(Amino_acids2, Amino_acids2))
# Amino_acids2_ = [i[0] + i[1] for i in Amino_acids2_]
Amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
               'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y']


def OE(seq):
    # 把氨基酸映射到整数；不在表里的统一用 0
    code_map = {aa: idx + 1 for idx, aa in enumerate(Amino_acids)}
    return [code_map.get(aa.upper(), 0) for aa in seq]


def AAC(seq):
    vec = np.zeros(21, dtype=np.float32)
    seq_len = len(seq)
    for aa in seq.upper():
        if aa in Amino_acids:
            vec[Amino_acids.index(aa)] += 1
    return vec / seq_len
    "AAC：21D"


def GGAP(seqs):
    GGAP_feature = []
    num = 0
    for i in range(len(seq) - 3):
        GGAP_feature.append((seq[i] + seq[i + 3]))

    seqs_.append([GGAP_feature.count(i) / (len(seq) - 3) for i in Amino_acids_])

    return seqs_


'上面是OE编码的结果：100D'


def DPC(seq, include_X=False):
    """
    顺序有关的二肽组成 (400 维)。
    include_X: 是否把 'X' 算进 21 类氨基酸，默认 False（20 类）
    """
    AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
           'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'] + (['X'] if include_X else [])
    n = len(AAs)
    vec = np.zeros(n * n, dtype=np.float32)

    # 二肽到索引映射
    dp2idx = {aa1 + aa2: i for i, (aa1, aa2) in enumerate(product(AAs, repeat=2))}

    seq = seq.upper()
    denom = max(len(seq) - 1, 1)
    for i in range(len(seq) - 1):
        key = seq[i:i + 2]
        if key in dp2idx:
            vec[dp2idx[key]] += 1
    return vec / denom


'DPC:231D'


def AAE_1(fastas):
    length = float(len(fastas))
    amino_acids = dict.fromkeys(letters, 0)
    encodings = []
    for AA in amino_acids:
        hits = [a.start() for a in list(re.finditer(AA, fastas))]
        p_prev = 0
        p_next = 1
        sum = 0
        while p_next < len(hits):
            distance = (hits[p_next] - hits[p_prev]) / length
            sum += distance * math.log(distance, 2)
            p_prev = p_next
            p_next += 1
        amino_acids[AA] = -sum
        encodings.append(amino_acids[AA])
    return encodings


#
def AAE(fastas):
    # encodings = []
    fastas_NT5 = "%s" % fastas[:5]
    fastas_CT5 = "%s" % fastas[-5:]
    encodings_full = AAE_1(fastas)
    encodings_CT5 = AAE_1(fastas_CT5)
    encodings_NT5 = AAE_1(fastas_NT5)
    # encodings.append(encodings_full + encodings_NT5 + encodings_CT5)
    return encodings_full + encodings_NT5 + encodings_CT5


# 没问题
def AAI(gene):
    with open("feature_file/AAI.txt") as f:
        records = f.readlines()[1:]
    AAI = []
    for i in records:
        array = i.rstrip().split()[1:] if i.rstrip() != '' else None
        AAI.append(array)
    AAI = np.array(
        [float(AAI[i][j]) for i in range(len(AAI)) for j in range(len(AAI[i]))]).reshape((14, 21))
    AAI = AAI.transpose()
    GENE_BE = {}
    AA = 'ACDEFGHIKLMNPQRSTWYV*'
    for i in range(len(AA)):
        GENE_BE[AA[i]] = i
    n = len(gene)
    gene_array = np.zeros((n, 14))
    for i in range(n):
        if gene[i] in GENE_BE:
            gene_array[i] = AAI[(GENE_BE[gene[i]])]
        else:
            gene_array[i] = np.zeros(14)
    return gene_array


def PC6_embedding(seq, max_len=100):
    f = open('feature_file/6-pc')
    text = f.read()
    f.close()
    text = text.split('\n')
    while '' in text:
        text.remove('')
    text = text[1:]
    AAI_dict = {}
    for each_line in text:
        temp = each_line.split(' ')
        while '' in temp:
            temp.remove('')
        for i in range(1, len(temp)):
            temp[i] = float(temp[i])
        AAI_dict[temp[0]] = temp[1:]
    AAI_dict['X'] = np.zeros(6)
    all_embeddings = []
    for each_seq in seq:
        # temp_embeddings = []
        # for each_char in each_seq:
        all_embeddings.append(AAI_dict[each_seq])
    # if max_len > len(each_seq):
    #     zero_padding = np.zeros((max_len - len(each_seq), 6))
    #     data_pad = np.vstack((temp_embeddings, zero_padding))
    # elif max_len == len(each_seq):
    #     data_pad = temp_embeddings
    # else:
    #     data_pad = temp_embeddings[:max_len]
    # all_embeddings.append(data_pad)
    all_embeddings = np.array(all_embeddings)
    return all_embeddings


def BLOSUM62(gene):
    with open("feature_file/blosum62.txt") as f:
        records = f.readlines()[1:]
    blosum62 = []
    for i in records:
        array = i.rstrip().split() if i.rstrip() != '' else None
        blosum62.append(array)
    blosum62 = np.array(
        [float(blosum62[i][j]) for i in range(len(blosum62)) for j in range(len(blosum62[i]))]).reshape((20, 21))
    blosum62 = blosum62.transpose()
    GENE_BE = {}
    AA = 'ARNDCQEGHILKMFPSTWYV'
    for i in range(len(AA)):
        GENE_BE[AA[i]] = i
    n = len(gene)
    gene_array = np.zeros((n, 20))
    for i in range(n):
        if gene[i] in GENE_BE:
            gene_array[i] = blosum62[(GENE_BE[gene[i]])]
        else:
            gene_array[i] = np.zeros(20)
    return gene_array


def BE(gene):
    with open("feature_file/BE.txt") as f:
        records = f.readlines()[1:]
    BE = []
    for i in records:
        array = i.rstrip().split() if i.rstrip() != '' else None
        BE.append(array)
    BE = np.array(
        [float(BE[i][j]) for i in range(len(BE)) for j in range(len(BE[i]))]).reshape((20, 20))
    BE = BE.transpose()
    AA = 'ACDEFGHIKLMNPQRSTWYV'
    GENE_BE = {}
    for i in range(len(AA)):
        GENE_BE[AA[i]] = i

    m = len(gene)
    # n = max(len(seq) for seq in train_sequences)
    gene_array = np.zeros((m, 20))
    for i in range(m):
        if gene[i] in GENE_BE:
            # 标准氨基酸：使用BE矩阵对应行
            gene_array[i] = BE[GENE_BE[gene[i]]]
        else:
            # 未知氨基酸（如X）：保留全零或其他策略
            gene_array[i] = np.zeros(20)

    return gene_array


GROUP = {
    'alphaticr': 'GAVLMI',
    'aromatic': 'FYW',
    'postive': 'KRH',
    'negative': 'DE',
    'neutral': 'STCPNQ'
}
GAP = 5  # 最大间隔
GROUP_KEYS = list(GROUP.keys())
GROUP_INDEX = {aa: k for k, v in GROUP.items() for aa in v}


def generateGroupPairs(groupKey):
    gPair = {}
    for key1 in groupKey:
        for key2 in groupKey:
            gPair[key1 + '.' + key2] = 0
    return gPair


def CKSAAGP(fastas, gap=5, **kw):
    group = {
        'alphaticr': 'GAVLMI',
        'aromatic': 'FYW',
        'postivecharger': 'KRH',
        'negativecharger': 'DE',
        'uncharger': 'STCPNQ'
    }

    AA = 'ARNDCQEGHILKMFPSTWYV'

    groupKey = group.keys()

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    gPairIndex = []
    for key1 in groupKey:
        for key2 in groupKey:
            gPairIndex.append(key1 + '.' + key2)

    encodings = []
    header = ['#']
    for g in range(gap + 1):
        for p in gPairIndex:
            header.append(p + '.gap' + str(g))
    encodings.append(header)

    sequence = re.sub('-', '', fastas)
    code = []
    for g in range(gap + 1):
        gPair = generateGroupPairs(groupKey)
        sum = 0
        for p1 in range(len(sequence)):
            p2 = p1 + g + 1
            if p2 < len(sequence) and sequence[p1] in AA and sequence[p2] in AA:
                gPair[index[sequence[p1]] + '.' + index[sequence[p2]]] = gPair[index[sequence[p1]] + '.' + index[
                    sequence[p2]]] + 1
                sum = sum + 1
        if sum == 0:
            for gp in gPairIndex:
                code.append(0)
        else:
            for gp in gPairIndex:
                code.append(gPair[gp] / sum)
    return code


'cks:149D'
GROUP = ['AGV', 'ILFP', 'YMTS', 'HNQW', 'RK', 'DE', 'C']
TRANS = {aa: str(i) for i, grp in enumerate(GROUP) for aa in grp}
BASE = 7
K = 3
N_TRIS = BASE ** K  # 343


def KMER(seq: str) -> np.ndarray:
    """
    输入：单条蛋白序列
    输出：343 维 7-进制 3-mer 频率向量
    """
    # 2. 转 7-进制字符串
    trans_seq = ''.join(TRANS.get(aa, '6') for aa in seq.upper())

    # 3. 统计 343 类 3-mer 出现次数
    cnt = np.zeros(N_TRIS, dtype=np.float32)
    for i in range(len(trans_seq) - K + 1):
        kmer = trans_seq[i:i + K]
        idx = int(kmer, BASE)  # 7 进制 → 0-3ES42
        cnt[idx] += 1

    # 4. 归一化到频率
    if len(trans_seq) >= K:
        cnt /= (len(trans_seq) - K + 1)
    return cnt  # (343,)


# 7 种物理化学属性分组
# GROUPS = [
#     {'hydrophobicity': 'RKEDQN'},  # low
#     {'hydrophobicity': 'GASTPHY'}, # medium
#     {'hydrophobicity': 'CLVIMFW'}, # high
#     {'volume':         'GASTPDC'},
#     {'volume':         'NVEQIL'},
#     {'volume':         'MHKFRYW'},
#     {'polarity':       'LIFWCMVY'},
#     {'polarity':       'PATGS'},
#     {'polarity':       'HQRKNED'},
#     {'polarizability': 'GASDT'},
#     {'polarizability': 'CPNVEQIL'},
#     {'polarizability': 'KMHFRYW'},
#     {'charge':         'KR'},
#     {'charge':         'ANCQGHILMFPSTWYV'},
#     {'charge':         'DE'},
#     {'ss':             'EALMQKRH'},
#     {'ss':             'VIYCWFT'},
#     {'ss':             'GNPSD'},
#     {'solvent':        'ALFCGIVW'},
#     {'solvent':        'RKQEND'},
#     {'solvent':        'MSPTHY'}
# ]

import numpy as np

GROUPS = {'alphaticr': 'GAVLMI', 'aromatic': 'FYW', 'postive': 'KRH', 'negative': 'DE', 'neutral': 'STCPNQ'}
GROUP_KEYS = list(GROUPS.keys())
GROUP_INDEX = {aa: k for k, v in enumerate(GROUPS.values()) for aa in v}


def _count(seq, group):
    return sum(aa in group for aa in seq) / len(seq)


def _transition(seq, group1, group2):
    pairs = [seq[i:i + 2] for i in range(len(seq) - 1)]
    cnt = sum((a in group1 and b in group2) or (a in group2 and b in group1) for a, b in pairs)
    return cnt / max(1, len(pairs))


def _distribution(seq, group):
    pos = [i + 1 for i, aa in enumerate(seq) if aa in group]
    if not pos: return [0] * 5
    quantile = [0.25, 0.5, 0.75, 1.0]
    return [(np.percentile(pos, q * 100) / len(seq)) for q in quantile]


def CTD(seq: str) -> np.ndarray:
    seq = seq.upper()
    code = []
    for g1, g2, g3 in zip(GROUPS[::3], GROUPS[1::3], GROUPS[2::3]):
        # Composition
        c1, c2, c3 = _count(seq, g1), _count(seq, g2), _count(seq, g3)
        code.extend([c1, c2, c3])
        # Transition
        t12, t13, t23 = _transition(seq, g1, g2), _transition(seq, g1, g3), _transition(seq, g2, g3)
        code.extend([t12, t13, t23])
        # Distribution
        code.extend(_distribution(seq, g1) + _distribution(seq, g2) + _distribution(seq, g3))
    return np.array(code)  # 147 维（可扩展为 343）


#
# def BLOSUM62_embedding(seq, max_len=100):
#     f = open('feature_file/blosum62.txt')
#     text = f.read()
#     f.close()
#     text = text.split('\n')
#     while '' in text:
#         text.remove('')
#     cha = text[0].split(' ')
#     while '' in cha:
#         cha.remove('')
#     index = []
#     for i in range(1, len(text)):
#         temp = text[i].split(' ')
#         while '' in temp:
#             temp.remove('')
#         for j in range(len(temp)):
#             temp[j] = float(temp[j])
#         index.append(temp)
#     index = np.array(index)
#     BLOSUM62_dict = {}
#     for j in range(len(cha)):
#         BLOSUM62_dict[cha[j]] = index[:, j]
#     all_embeddings = []
#     for each_seq in seq:
#         temp_embeddings = []
#         for each_char in each_seq:
#             temp_embeddings.append(BLOSUM62_dict[each_char])
#         if max_len > len(each_seq):
#             zero_padding = np.zeros((max_len - len(each_seq), 23))
#             data_pad = np.vstack((temp_embeddings, zero_padding))
#         elif max_len == len(each_seq):
#             data_pad = temp_embeddings
#         else:
#             data_pad = temp_embeddings[:max_len]
#         all_embeddings.append(data_pad)
#     all_embeddings = np.array(all_embeddings)
#     return torch.from_numpy(all_embeddings).float()
#
# if 'BLO2' in featurelist:
#     maxabs = preprocessing.MinMaxScaler()
#     n1 = len(train_sequences)
#     BLO2_TRAIN=BLOSUM62_embedding(train_sequences)
#     BLO2_TRAIN =BLO2_TRAIN.reshape(-1, 1)
#     BLO2_TRAIN = maxabs.fit_transform(BLO2_TRAIN)
#     BLO2_TRAIN = BLO2_TRAIN.reshape(n1, 100,-1)
#     BLO2_TEST =BLO2_TEST.reshape(-1, 1)
#     BLO2_TEST = maxabs.fit_transform(BLO2_TEST)
#     a_train_2d.append(BLO2_TRAIN)
#
#
# if 'AAI' in featurelist:
#     def AAI(gene):
#         with open("feature_file/AAI.txt") as f:
#             records = f.readlines()[1:]
#         AAI = []
#         for i in records:
#             array = i.rstrip().split()[1:] if i.rstrip() != '' else None
#             AAI.append(array)
#         AAI = np.array(
#             [float(AAI[i][j]) for i in range(len(AAI)) for j in range(len(AAI[i]))]).reshape((14, 21))
#         AAI = AAI.transpose()
#         GENE_BE = {}
#         AA = 'ACDEFGHIKLMNPQRSTWYV*'
#         for i in range(len(AA)):
#             GENE_BE[AA[i]] = i
#         n = len(gene)
#         gene_array = np.zeros((100, 14))
#         for i in range(n):
#             gene_array[i] = AAI[(GENE_BE[gene[i]])]
#         return gene_array
#
#
#     n = len(train_sequences)
#     maxabs = preprocessing.MinMaxScaler()
#     x_AAI = np.zeros((n, 100, 14))  # 3880*100*20
#     for i in range(n):
#         x_AAI[i] = AAI(train_sequences[i])
#     # x_AAI = x_AAI.reshape(-1, 1)
#     # x_AAI = maxabs.fit_transform(x_AAI)
#     # x_AAI = x_AAI.reshape(n, 100*14)
#     # x_AAI2 = x_AAI2.reshape(-1, 1)
#     # x_AAI2 = maxabs.fit_transform(x_AAI2)
#     # x_AAI2 = x_AAI2.reshape(n2, 100*14)
#     a_train_2d.append(x_AAI)
#

"""PAAC"""
Amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
               'R', 'S', 'T', 'V', 'W', 'Y']


# 40D
def PAAC(seq, max_len=100):
    f = open('feature_file/PAAC.txt')
    text = f.read()
    f.close()
    text = text.split('\n')
    while '' in text:
        text.remove('')
    cha = text[0].split('\t')
    while '' in cha:
        cha.remove('')
    cha = cha[1:]
    index = []
    for i in range(1, len(text)):
        temp = text[i].split('\t')
        while '' in temp:
            temp.remove('')
        temp = temp[1:]
        for j in range(len(temp)):
            temp[j] = float(temp[j])
        index.append(temp)
    index = np.array(index)
    AAI_dict = {}
    for j in range(len(cha)):
        AAI_dict[cha[j]] = index[:, j]
    AAI_dict['X'] = np.zeros(3)
    temp_embeddings = []
    for each_char in seq:
        temp_embeddings.append(AAI_dict[each_char])
    if max_len > len(seq):
        zero_padding = np.zeros((max_len - len(seq), 3))
        data_pad = np.vstack((temp_embeddings, zero_padding))
    elif max_len == len(seq):
        data_pad = temp_embeddings
    else:
        data_pad = temp_embeddings[:max_len]
    # all_embeddings = np.array(all_embeddings)
    return data_pad


# def PSAAC(seqs):
#     seqs_ = []
#     PSAAC_profile_forward = []
#     PSAAC_profile_backward = []
#     forward_seq = []
#     backward_seq = []
#     i = 0
#     for seq in seqs:
#         forward_seq.append(list(seq[:5]))
#         backward_seq.append(list(seq[-5:]))
#
#     for position in range(5):
#         PSAAC_profile_forward.append(
#             [list(np.array(forward_seq)[:, position]).count(amino) / len(seqs) for amino in Amino_acids])
#
#     for position in range(5):
#         PSAAC_profile_backward.append(
#             [list(np.array(backward_seq)[:, position]).count(amino) / len(seqs) for amino in Amino_acids])
#
#     for seq in forward_seq:
#         num = 0
#         new_seq = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#         for amino in seq:
#             index_ = Amino_acids.index(amino)
#             new_seq[index_] = np.array(PSAAC_profile_forward)[num, index_]
#             num += 1
#
#         seqs_.append(new_seq)
#
#     for seq in backward_seq:
#         num = 0
#         new_seq = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#         for amino in seq:
#             index_ = Amino_acids.index(amino)
#             new_seq[index_] = np.array(PSAAC_profile_backward)[num, index_]
#             num += 1
#
#         seqs_[i].extend(new_seq)
#         i += 1
#     return seqs_
# if 'PSAAC' in featurelist:
#     psaac_train = np.array(PSAAC(train_sequences))
#     a_train.append(psaac_train)
# def PSAAC2(seqs):
#     seqs_ = []
#     PSAAC_profile_forward = []
#     PSAAC_profile_backward = []
#     forward_seq = []
#     backward_seq = []
#     i = 0
#     for seq in seqs:
#         forward_seq.append(list(seq[:5]))
#         backward_seq.append(list(seq[-5:]))
#
#     for position in range(5):
#         PSAAC_profile_forward.append(
#             [list(np.array(forward_seq)[:, position]).count(amino) / len(seqs) for amino in Amino_acids2])
#
#     for position in range(5):
#         PSAAC_profile_backward.append(
#             [list(np.array(backward_seq)[:, position]).count(amino) / len(seqs) for amino in Amino_acids2])
#
#     for seq in forward_seq:
#         num = 0
#         new_seq = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#         for amino in seq:
#             index_ = Amino_acids2.index(amino)
#             new_seq[index_] = np.array(PSAAC_profile_forward)[num, index_]
#             num += 1
#
#         seqs_.append(new_seq)
#
#     for seq in backward_seq:
#         num = 0
#         new_seq = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#         for amino in seq:
#             index_ = Amino_acids2.index(amino)
#             new_seq[index_] = np.array(PSAAC_profile_backward)[num, index_]
#             num += 1
#
#         seqs_[i].extend(new_seq)
#         i += 1
#     return seqs_
# if 'PSAAC2' in featurelist:
#     psaac_train = np.array(PSAAC2(train_sequences))
#     a_train.append(psaac_train)
#
# def AAI_embedding(seq, max_len=100):
#     f = open('feature_file/AAindex_.txt')
#     text = f.read()
#     f.close()
#     text = text.split('\n')
#     while '' in text:
#         text.remove('')
#     cha = text[0].split('\t')
#     while '' in cha:
#         cha.remove('')
#     cha = cha[1:]
#     index = []
#     for i in range(1, len(text)):
#         temp = text[i].split('\t')
#         while '' in temp:
#             temp.remove('')
#         temp = temp[1:]
#         for j in range(len(temp)):
#             temp[j] = float(temp[j])
#         index.append(temp)
#     index = np.array(index)
#     AAI_dict = {}
#     for j in range(len(cha)):
#         AAI_dict[cha[j]] = index[:, j]
#     AAI_dict['X'] = np.zeros(531)
#     all_embeddings = []
#     for each_seq in seq:
#         temp_embeddings = []
#         for each_char in each_seq:
#             temp_embeddings.append(AAI_dict[each_char])
#         if max_len > len(each_seq):
#             zero_padding = np.zeros((max_len - len(each_seq), 531))
#             data_pad = np.vstack((temp_embeddings, zero_padding))
#         elif max_len == len(each_seq):
#             data_pad = temp_embeddings
#         else:
#             data_pad = temp_embeddings[:max_len]
#         all_embeddings.append(data_pad)
#     all_embeddings = np.array(all_embeddings)
#     return torch.from_numpy(all_embeddings).float()
#
# if 'AAI2' in featurelist:
#     maxabs = preprocessing.MinMaxScaler()
#     n1 = len(train_sequences)
#     AAI_TRAIN=AAI_embedding(train_sequences)
#     AAI_TRAIN =AAI_TRAIN.reshape(-1, 1)
#     AAI_TRAIN = maxabs.fit_transform(AAI_TRAIN)
#     AAI_TRAIN = AAI_TRAIN.reshape(n1, 100, -1)
#     AAI_TEST = maxabs.fit_transform(AAI_TEST)
#     a_train_2d.append(AAI_TRAIN)
#
#
#
# # if 'CTD' in featurelist:
# #     def CTD(pep):  # Chain-Transition-Ditribution feature
# #         feature = []
# #         name = []
# #         for seq in pep:
# #             protein.ReadProteinSequence(seq)
# #             ctd = protein.GetCTD()
# #             feature.append(list(ctd.values()))
# #             name = list(ctd.keys())
# #         return feature, name
# def AAT(sequence, aatdic=readAAT("feature_file/aat_minmaxscaler_general.txt"),
#         avg=1):  # return AAT features for the peptides
#
#     feature = []
#     sequence = sequence.upper()  # 转换为大写
#     score = []
#     count = 0
#
#     for i in range(len(sequence) - 2):  # 遍历序列生成所有可能的三联体
#         try:
#             # 尝试从字典中获取三联体的分数
#             score.append(round(float(aatdic[sequence[i:i + 3]]), 4))
#             count += 1
#         except KeyError:  # 如果三联体不在字典中，则添加 -1 分数
#             score.append(float(-1))
#             count += 1
#     score = np.array(score, dtype=np.float32)
#     # 如果avg为0，则直接返回分数列表；如果avg为1，则计算平均分数
#     if int(avg) == 0:
#         feature.extend(score)
#     elif int(avg) == 1:
#         average_score = np.mean(score) / count if count != 0 else 0
#         feature.append(round(float(average_score), 4))
#     return feature
#
#     # feature = []
#     # for a in pep:
#     #     if int(avg) == 0:
#     #         # print(a)
#     #         score = []
#     #         count = 0
#     #         for i in range(0, len(a) - 2):
#     #             try:
#     #                 score.append(round(float(aatdic[a[i:i + 3]]), 4))
#     #                 # score += float(aapdic[a[i:i + 3]])
#     #                 count += 1
#     #             except KeyError:
#     #                 # print(a[i:i + 3])
#     #                 score.append(float(-1))
#     #                 # score += -1
#     #                 count += 1
#     #                 continue
#     #         # averagescore = score / count
#     #         feature.append(score)
#     #     if int(avg) == 1:
#     #         score = 0
#     #         count = 0
#     #         for i in range(0, len(a) - 2):
#     #             try:
#     #                 score += float(aatdic[a[i:i + 3]])
#     #                 count += 1
#     #             except KeyError:
#     #                 score += -1
#     #                 count += 1
#     #                 continue
#     #         # print(a, score)
#     #         if count != 0:
#     #             averagescore = score / count
#     #         else:
#     #             averagescore = 0
#     #         feature.append(round(float(averagescore), 4))
#     # return feature


#
# def AAP(sequence, aapdic=readAAP("feature_file/aap_minmaxscaler_general.txt"),
#         avg=1):  # return AAP features for the peptides
#     feature = []
#     sequence = sequence.upper()  # 转换为大写
#     score = []
#     count = 0
#
#     for i in range(len(sequence) - 1):  # 遍历序列生成所有可能的二联体
#         try:
#             # 尝试从字典中获取二联体的分数
#             score.append(round(float(aapdic[sequence[i:i + 2]]), 4))
#             count += 1
#         except KeyError:  # 如果二联体不在字典中，则添加 -1 分数
#             score.append(float(-1))
#             count += 1
#     score = np.array(score, dtype=np.float32)
#     # 如果avg为0，则直接返回分数列表；如果avg为1，则计算平均分数
#     if int(avg) == 0:
#         feature.extend(score)
#     elif int(avg) == 1:
#         average_score = np.mean(score) / count if count != 0 else 0
#         feature.append(round(float(average_score), 4))
#     return feature


# def PC6_embedding(seq, max_len=100):
#     f = open('2d_features/6-pc')
#     text = f.read()
#     f.close()
#     text = text.split('\n')
#     while '' in text:
#         text.remove('')
#     text = text[1:]
#     AAI_dict = {}
#     for each_line in text:
#         temp = each_line.split(' ')
#         while '' in temp:
#             temp.remove('')
#         for i in range(1, len(temp)):
#             temp[i] = float(temp[i])
#         AAI_dict[temp[0]] = temp[1:]
#     AAI_dict['X'] = np.zeros(6)
#     all_embeddings = []
#     for each_seq in seq:
#         temp_embeddings = []
#         for each_char in each_seq:
#             temp_embeddings.append(AAI_dict[each_char])
#         if max_len > len(each_seq):
#             zero_padding = np.zeros((max_len - len(each_seq), 6))
#             data_pad = np.vstack((temp_embeddings, zero_padding))
#         elif max_len == len(each_seq):
#             data_pad = temp_embeddings
#         else:
#             data_pad = temp_embeddings[:max_len]
#         all_embeddings.append(data_pad)
#     all_embeddings = np.array(all_embeddings)
#     return torch.from_numpy(all_embeddings).float()
# if 'PC6' in featurelist:
#     maxabs = preprocessing.MinMaxScaler()
#     n1 = len(train_sequences)
#     PC6_TRAIN=PC6_embedding(train_sequences)
#     # PC6_TRAIN =PC6_TRAIN.reshape(-1, 1)
#     # PC6_TRAIN = maxabs.fit_transform(PC6_TRAIN)
#     # PC6_TRAIN = PC6_TRAIN.reshape(n1, 100 * 6)
#     # PC6_TEST =PC6_TEST.reshape(-1, 1)
#     # PC6_TEST = maxabs.fit_transform(PC6_TEST)
#     # PC6_TEST = PC6_TEST.reshape(n2, 100 * 6)
#     a_train_2d.append(PC6_TRAIN)
#
# def Skip(seq, skip):
#     element = []
#     for i in range(len(seq) - skip - 1):
#         element.append(seq[i] + seq[i + skip + 1])
#     return element
#
# def ASDC(seqs):
#     seqs_ = []
#     for seq in seqs:
#         ASDC_feature = []
#         skip = 0
#         for i in range(len(seq)):
#             ASDC_feature.extend(Skip(seq, skip))
#             skip += 1
#         seqs_.append([ASDC_feature.count(i) / len(ASDC_feature) for i in Amino_acids_])
#     return seqs_
#
# if 'asdc' in featurelist:
#     asdc_train = np.array(ASDC(train_sequences))
#     a_train.append(asdc_train)
def ZS(seqs):
    zscale = {
        'A': [0.24, -2.32, 0.60, -0.14, 1.30],  # A
        'C': [0.84, -1.67, 3.71, 0.18, -2.65],  # C
        'D': [3.98, 0.93, 1.93, -2.46, 0.75],  # D
        'E': [3.11, 0.26, -0.11, -0.34, -0.25],  # E
        'F': [-4.22, 1.94, 1.06, 0.54, -0.62],  # F
        'G': [2.05, -4.06, 0.36, -0.82, -0.38],  # G
        'H': [2.47, 1.95, 0.26, 3.90, 0.09],  # H
        'I': [-3.89, -1.73, -1.71, -0.84, 0.26],  # I
        'K': [2.29, 0.89, -2.49, 1.49, 0.31],  # K
        'L': [-4.28, -1.30, -1.49, -0.72, 0.84],  # L
        'M': [-2.85, -0.22, 0.47, 1.94, -0.98],  # M
        'N': [3.05, 1.62, 1.04, -1.15, 1.61],  # N
        'P': [-1.66, 0.27, 1.84, 0.70, 2.00],  # P
        'Q': [1.75, 0.50, -1.44, -1.34, 0.66],  # Q
        'R': [3.52, 2.50, -3.50, 1.99, -0.17],  # R
        'S': [2.39, -1.07, 1.15, -1.39, 0.67],  # S
        'T': [0.75, -2.18, -1.12, -1.46, -0.40],  # T
        'V': [-2.59, -2.64, -1.54, -0.85, -0.02],  # V
        'W': [-4.36, 3.94, 0.59, 3.44, -1.59],  # W
        'Y': [-2.54, 2.44, 0.43, 0.04, -1.47],  # Y
        '-': [0.00, 0.00, 0.00, 0.00, 0.00],  # -
    }
    # encodings = []
    encoding = []
    for aa in seqs:
        encoding.append(zscale[aa])

    # encodings.append(encoding)

    return encoding


AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
      'R', 'S', 'T', 'V', 'W', 'Y', 'X']
CODON_WT = {'A': 4, 'C': 2, 'D': 2, 'E': 2, 'F': 2, 'G': 4, 'H': 2, 'I': 3,
            'K': 2, 'L': 6, 'M': 1, 'N': 2, 'P': 4, 'Q': 2, 'R': 6, 'S': 6,
            'T': 4, 'V': 4, 'W': 1, 'Y': 2, 'X': 1}


def DDE(seq: str) -> np.ndarray:
    seq = re.sub('-', '', seq.upper())
    L = max(1, len(seq) - 1)  # 避免除零
    # 1. 二肽计数
    # cnt = np.zeros(400, dtype=float)
    cnt = np.zeros(441, dtype=float)
    for i in range(len(seq) - 1):
        try:
            idx = AA.index(seq[i]) * 20 + AA.index(seq[i + 1])
            cnt[idx] += 1
        except ValueError:  # 处理未知氨基酸X
            idx = 0  # 将X映射到索引0
            cnt[idx] += 1
        # except ValueError:  # 处理未知氨基酸X
        #     idx = AA.index('X') * 20 + AA.index('X')  # 将X映射到自身
        #     cnt[idx] += 1
    cnt /= cnt.sum() if cnt.sum() else 1
    # 2. 期望 & 方差
    exp = np.array([(CODON_WT[a] / 61) * (CODON_WT[b] / 61) for a in AA for b in AA])
    var = exp * (1 - exp) / L
    var[var == 0] = 1e-12  # 防除零
    # 3. 标准化偏差
    dde = (cnt - exp) / np.sqrt(var)
    return dde  # (400,)


# if len(a_train) > 1:
#     combined_train = np.concatenate(a_train, axis=1)
# else:
#     combined_train = a_train[0]
#
# # QSOrder = iFeatureOmegaCLI.iProtein(train_sequences)
# # QSOrder.get_descriptor("QSOrder")
# #
# # GTPC = iFeatureOmegaCLI.iProtein(test_sequences)
# # GTPC.get_descriptor("GTPC type 2")
# # GTPC.encodings = GTPC.encodings.reset_index(drop=True)
# # QSOrder.encodings = QSOrder.encodings.reset_index(drop=True)
#
#
#
# if len(a_train_2d) > 0:
#     train_2d = np.concatenate(a_train_2d, axis=2)
# else:
#     train_2d=None


# ---------- 编码器表 ----------
ENCODERS = {
    # 'OE'  : OE,
    'AAC': AAC,
    'DPC': DPC,
    'AAE': AAE,
    'BE': BE,
    'Blosum62': BLOSUM62,
    'AAI': AAI,
    'PC6': PC6_embedding
    # 以后直接往这里加
}

def read_and_extract(filename):
    sequences = []  # 存储序列的列表

    with open(filename, 'r') as file:
        sequence = ''
        for line in file:
            if line.startswith('>'):  # 假设标签以 '>' 开头
                pass
            else:
                sequences.append(line.strip())
    return sequences

def combinefeature(input_file, encoders, pad_len=None):
    # protein_list = 'DATA/protein sequences.fasta'
    '所有的序列为train_sequence,test_sequence'
    peptide_sequences = []  # 3880
    # test_sequences=[]#970
    sequences = read_and_extract(input_file)
    peptide_sequences.extend(sequences)
    MAX_LENGTH = max(
        len(str(seq)) for seq in peptide_sequences)  # 强制转为字符串类型避免报错‌:ml-citation{ref="1" data="citationList"}
    print(f"第二列最长序列长度为：{MAX_LENGTH}")

    lengths = np.array([len(str(seq)) for seq in peptide_sequences])

    # 定义区间
    bins = [0, 50, 100, 150, 200, 300, 400, 500, 1000,2000,3000,4000]
    labels = [f"{bins[i]}-{bins[i + 1]}" for i in range(len(bins) - 1)]

    # 统计
    hist, _ = np.histogram(lengths, bins=bins)
    for rng, cnt in zip(labels, hist):
        print(f"{rng}: {cnt} 条")
    # records = {}
    feat2_dict, feat3_dict = {}, {}
    feature_dimensions_1d = {}
    feature_dimensions_2d = {}
    for record in SeqIO.parse(input_file, "fasta"):
        pid = record.id
        seq = str(record.seq).upper()
        # seq = 'XXXACDEF'
        try:
            feat2_list, feat3_list = [], []# list[int]  → 长度 = 序列长度
            for name, func in encoders.items():
                feat = np.asarray(func(seq), dtype=np.float32)



                # 统一长度（可选）
                if pad_len is not None:
                    if feat.ndim == 1:
                        feat = np.pad(feat, (0, max(0, pad_len - len(feat))))[:pad_len]
                    else:
                        # 二维/三维自行处理
                        pass
                if feat.ndim == 1:
                    feature_dimensions_1d[name] = feat.shape
                    # feat = np.pad(feat, (0, max(0, 100 - len(feat))), 'constant')[:100]
                    feat2_list.append(feat)
                elif feat.ndim == 2:
                    # feat = np.pad(feat, (0, max(0, 100 - len(feat))), 'constant')[:100]
                    feat = np.pad(feat,
                                  ((0, max(0, 100 - feat.shape[0])),  # 行维补 0
                                   (0, 0)),  # 列维不补
                                  'constant')
                    feature_dimensions_2d[name] = feat.shape
                    feat3_list.append(feat)

                    # 拼接
            if feat2_list:
                    feat2_dict[pid] = np.concatenate(feat2_list, axis=-1)  # (L, C2)
            if feat3_list:
                    feat3_dict[pid] = np.concatenate(feat3_list, axis=-1)  # (L, C3)
        except Exception as e:
            print(f"Skip {pid}: {e}")
            import traceback
            traceback.print_exc()

    df_dimensions = pd.DataFrame(columns=['Feature Name', 'Dimensions'])
    for name, shape in feature_dimensions_1d.items():
        df_dimensions = df_dimensions._append({'Feature Name': name, 'Dimensions': str(shape)}, ignore_index=True)
    for name, shape in feature_dimensions_2d.items():
        df_dimensions = df_dimensions._append({'Feature Name': name, 'Dimensions': str(shape)}, ignore_index=True)
    # table_data = []
    # # 填充列表
    # for name, shape in df_dimensions.items():
    #     table_data.append([name, str(shape)])

    # 使用 tabulate 库创建表格
    # table = tabulate(table_data, headers=["Feature Name", "Dimensions"], tablefmt="fancyGrid")
    # print(table)

    #print(f"Total feature dimensions: {total_dimensions}")
    total_dimensions = 0
    # for name, func in encoders.items():
    #     feat = func("MKFLVLFLVLFLCLLGVGLQ")  # 示例序列
    #     feature_dimensions[name] = feat.size
    #     total_dimensions += np.prod(feat.shape)

    # df_dimensions = pd.DataFrame(columns=['Feature Name', 'Dimensions'])
    #
    # # 填充 DataFrame
    # for name, shape in feature_dimensions.items():
    #     df_dimensions = df_dimensions.append({'Feature Name': name, 'Dimensions': str(shape)}, ignore_index=True)
    #
    print(df_dimensions)
    return feat2_dict, feat3_dict

# 1. 切到数据目录
# os.chdir('/home/researchlab/Downloads/CODES/Data')
# print('PWD:', os.getcwd())
pept_emb_file = '../Feature Development/Input_dataset.fasta'
pept_hand_2dict,pept_hand_3dict = combinefeature(pept_emb_file, encoders={'AAC': AAC,
                                                                    'DPC': DPC, 'AAE':AAE, 'BE' :BE,
    'Blosum62' : BLOSUM62,
    'AAI' : AAI,
    'PC6' : PC6_embedding,
    'DDE' : DDE,
    # 'AAT' : AAT,
    # 'AAP' : AAP,
    #'PAAC' :PAAC,
    'KMER': KMER,
    #'CKSGGAP':CKSAAGP,
    #'ZS' : ZS,
    #'CTD' : CTD,
                                                       })
    # 以后直接往这里加)
all_values = list(pept_hand_2dict.values())   # 取全部value
# import numpy as np
X_union_2 = np.stack(all_values)


# 2. 读取主数据集
df_union = pd.read_csv('9_dataset_subset_unique.csv')
df_union_columns = df_union.columns

# 3. 构造 X, Y
# X_union = df_union.drop(['ID', 'Label_', 'Sequence'], axis=1)

X_union = np.load(r'D:\Desktop\anti-ahv\Pred-AHCP-main\Pred-AHCP-main\Feature Development\ahcp_ESM2_35m_2d.npy', allow_pickle=True).item()
# X_union = np.load(r'D:\Desktop\anti-ahv\Pred-AHCP-main\Pred-AHCP-main\Feature Development\ahcp_ESM2_650m_2d.npy', allow_pickle=True).item()

Y_union = np.ravel(df_union['Label'])
all_values = list(X_union.values())   # 取全部value
# import numpy as np
X_union= np.stack(all_values)
X_ALL = np.concatenate((X_union_2,X_union),axis=1)


from sklearn.preprocessing import MinMaxScaler
scaler1 = MinMaxScaler()
X_union = scaler1.fit_transform(X_ALL)
#653*480        # shape: (N_samples, L, feat_dim)
# 4. 归一化>1的列
# over1 = X_union.columns[X_union.max() > 1.0]
# X_union[over1] = X_union[over1] / 100

# 5. 训练/测试划分
from sklearn.model_selection import train_test_split
X_union_train, X_union_test, Y_union_train, Y_union_test = \
    train_test_split(X_union, Y_union, test_size=0.2, random_state=99)
np.save('X_origin_space.npy', X_union_test)
# 1. 对 DataFrame 做完全相同的划分（random_state 保持一致）
df_train, df_test = train_test_split(
    df_union,                # 整表拆分
    test_size=0.2,
    random_state=99,
    stratify=df_union['Label']
)
test_sequences = df_test['Sequence'].values
test_ids       = [f"seq_{idx}" for idx in df_test.index]   # 用原索引当 ID

# 3. 写入 FASTA
import time
timestamp = time.strftime("%Y%m%d_%H%M%S")
fasta_path = fr"D:\Desktop\anti-ahv\Pred-AHCP-main\Pred-AHCP-main\test_sequences_{timestamp}.fasta"

with open(fasta_path, 'w') as f:
    for seq_id, seq, label in zip(test_ids, test_sequences, df_test['Label']):
        f.write(f">{seq_id}|label={label}\n{seq}\n")
# with open(fasta_path, 'w') as f:
#     for seq_id, seq in zip(test_ids, test_sequences):
#         f.write(f">{seq_id}\n{seq}\n")

print(f"测试集 FASTA 已保存至：{fasta_path}")

import lightgbm as lgb
from sklearn.feature_selection import SelectFromModel
# ---------- 1. 用 LightGBM 快速选特征 ----------
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
# 基模型：轻量级、速度快，RFE 迭代不耗时
base_clf = LogisticRegression(max_iter=1000, solver='lbfgs')

# 保留前 30% 特征（可调）
n_feat = int(0.3 * X_union_train.shape[1])   # 例如 1745 -> 523
rfe = RFE(estimator=base_clf, n_features_to_select=n_feat, step=0.1)
rfe.fit(X_union_train, Y_union_train)
X_union_train = rfe.transform(X_union_train)
X_union_test = rfe.transform(X_union_test)      # 测试集禁止再 fit
X_union= rfe.transform(X_union)           # 交叉验证用
np.save('X_selected_space.npy', X_union_test)
# # 统一转换
# X_train_sel = selector.transform(X_union_train)
# X_test_sel  = selector.transform(X_union_test)
# X_union_sel = selector.transform(X_union)      # 用于交叉验证

# print(f'原始特征数: {X_union_train.shape[1]} -> 筛选后: {X_train_sel.shape[1]}')


# 6. 读取独立验证集
df_val = pd.read_csv('Validation.csv')[df_union_columns]
X_val = df_val.drop(['ID', 'Label', 'Sequence'], axis=1)
Y_val = np.ravel(df_val['Label'])

# 7. 定义分类器
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from xgboost import XGBClassifier
from catboost import CatBoostClassifier
# 'Random Forest': RandomForestClassifier(n_estimators=135, max_depth=15, min_samples_leaf=2, max_features=int(1265*0.5), random_state=42),
# 'Logistic Regression': LogisticRegression(max_iter=100),#无随机
# 'Support Vector Machine poly': SVC(kernel='poly', C=1.0, gamma='scale', random_state=42, probability=True),
# max_feat_range = range(50, 1266, 50)
# for mf in max_feat_range:


classifiers_5 = {

    'Random Forest':RandomForestClassifier(n_estimators=135, max_depth=14, min_samples_leaf=2,max_features=15, random_state=42),#93.1
    'Logistic Regression': LogisticRegression(max_iter=100),#无随机#92.4
    'Support Vector Machine poly': SVC(kernel='poly', C=0.4, gamma='scale', random_state=42, probability=True),#93.9
    'LGBM': lgb.LGBMClassifier(
            n_estimators=200,
             num_leaves=7,          # 调小叶子数
             min_data_in_leaf=5,
             max_depth =  6,
        # 每个叶子最少样本
             feature_fraction= 0.2,
             bagging_fraction= 0.2,
             learning_rate= 0.1,
             verbose=-1,             # 关键：关闭 LightGBM 日志
             random_state=42,
         ),#90.8

    'CatBoost':CatBoostClassifier(
        iterations=300,  # 迭代次数，可以根据数据量和训练时间进行调整
        learning_rate=0.14,  # 学习率，控制每次迭代的步长
        depth=6,  # 树的深度，控制模型的复杂度
        loss_function='Logloss',  # 损失函数，对于分类任务通常选择Logloss
        # custom_metric=['AUC'],
        l2_leaf_reg= 3,
        bagging_temperature=1,
        random_seed=42  # 随机种子，用于复现结果
    )#91.6
}

# 8. 训练+评估
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn import metrics
from sklearn.ensemble import VotingClassifier


summary = []
results = []
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
five = {'SVM': classifiers_5['Support Vector Machine poly'],
        'RF': classifiers_5['Random Forest'],
        'LR': classifiers_5['Logistic Regression'],
        'CatBoost': classifiers_5['CatBoost'],
        # 'LGBM':classifiers_5['LGBM']
        }
#
# # 2. 软投票集成
soft5 = VotingClassifier(estimators=list(five.items()),
                         voting='hard',
                         n_jobs=-1)
#
# ================= SHAP Analysis =================
"***************************************"
import shap
import matplotlib.pyplot as plt

rf_model = classifiers_5['Random Forest']
rf_model.fit(X_union_train, Y_union_train)

explainer = shap.TreeExplainer(rf_model)
shap_values = explainer.shap_values(X_union_train)
print("Train shape:", X_union_train.shape)
print("SHAP shape:", shap_values[0].shape)
# Summary plot
plt.figure()
shap.summary_plot(shap_values[:,:,1], X_union_train, show=False)
# shap.summary_plot(shap_values[1], X_union_train, max_display=20)

plt.tight_layout()
plt.savefig("shap_summary.png", dpi=600)
plt.close()

# Bar plot
plt.figure()
shap.summary_plot(shap_values[:, :, 1], X_union_train, plot_type="bar", show=False)
plt.tight_layout()
plt.savefig("shap_bar.png", dpi=600)
plt.close()
"***************************************"
#
# # # 3. 训练 + 评估（复用你原来的逻辑）
# soft5.fit(X_union_train, Y_union_train)
# y_pred = soft5.predict(X_union_test)
# # prob   = soft5.predict_proba(X_union_test)[:, 1]
# cv_scores = cross_val_score(
#     soft5,
#     X_union,
#     Y_union,
#     cv=cv,
#     scoring='accuracy',
#     n_jobs=-1
# )
# # # 4. 计算指标
# results.append({
#     'Classifier': 'Soft-Vote-5',
#     'Testing Accuracy': metrics.accuracy_score(Y_union_test, y_pred),
#     'Testing Precision': metrics.precision_score(Y_union_test, y_pred),
#     'Testing Recall': metrics.recall_score(Y_union_test, y_pred),
#     'Testing Specificity': metrics.recall_score(Y_union_test, y_pred, pos_label=0),
#     # 'Testing AUROC': metrics.roc_auc_score(Y_union_test, prob),
#     'Testing MCC': metrics.matthews_corrcoef(Y_union_test, y_pred),
#     'Testing F1 Score': metrics.f1_score(Y_union_test, y_pred),
#     'Mean Stratified CV Accuracy': cv_scores.mean(),
#     'Std Stratified CV Accuracy': cv_scores.std()
# })
#
# pd.set_option('display.max_columns', None)
# pd.set_option('display.width', None)
# # 5. 保存最终表格（含软投票一行）
# results_df = pd.DataFrame(results).round(3)
# # print(results_df)
# # import joblib
# #
# # joblib.dump(soft5, 'final_voting_model.pkl')
# #
# print(results_df.to_string(index=False))
# print('All done! results_df.csv & random_forest.pkl have been saved.')

# 计算shap值
# shap_values = explainer.shap_values(X_union_train)
# import joblib
#
# model = joblib.load('final_voting_model.pkl')
# probs = []
#
# for name, clf in model.named_estimators_.items():
#     p = clf.predict_proba(X_union_test)[:,1]
#     probs.append(p)
#
# import numpy as np
# X_model_space = np.vstack(probs).T
#
# X_CLASSIFIER = X_model_space.toarray() if hasattr(X_model_space, "toarray") else np.array(X_model_space)
# y_labels = np.array(Y_union_test) # 假设你的标签变量名是 y
#
# # 3. 保存为 .npy 文件
#
# np.save('X_model_features.npy',X_CLASSIFIER)
# np.save('y_labels.npy', y_labels)
#
# print("所有数据已成功保存为 .npy 格式！")