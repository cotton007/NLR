import argparse

def nlr_summary2dict(file):
    """
    将nlr_summary.csv中的候选基因及其信息存储到字典中
    :param file: nlr_summary.csv路径
    :return: {"基因名":"基因名\t evidence\t structure"}
    """
    result_dict = {}
    with open(file,'r') as f:
        line = True
        while line:
            line = f.readline()
            if not line:
                break
            line_list = line.strip().split()
            if line_list[0] != "gene":
                result_dict[line_list[0]] = line
    return result_dict


def seqlen(file):
    """
    输入fasta文件，统计每条序列的长度，存储为字典
    :param file:fasta文件路径
    :return:{"序列名"：序列长度}
    """
    result_dict = {}
    with open(file,'r') as f:
        line = True
        length = 0
        while line:
            line = f.readline()
            if not line:
                break
            line_list = line.strip().split()
            if line_list[0].startswith(">"):
                id = line_list[0][1:]
                length = 0
            else:
                length += len(line)
                result_dict[id] = length
    return result_dict


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--fasta",help='输入蛋白质fasta文件')
    parser.add_argument("-c", "--csv", help='输入nlr_summary.csv')
    parser.add_argument("-o", "--out", help='输出文件路径')
    args = parser.parse_args()

    peplen_dict = seqlen(args.fasta)
    # peplen_dict记录了每个蛋白质序列的长度:{"序列名"：序列长度}
    nlr_dict = nlr_summary2dict(args.csv)
    # nlr_dict记录了经鉴定的nlr候选基因，以及evidence和structure:{"基因名":"基因名\t evidence\t structure"}

    with open(args.out,'w') as f:
        f.write("gene\tevidence\tstructure\n")
        for i in nlr_dict:
            if peplen_dict[i] >= 150:
                f.write(nlr_dict[i])