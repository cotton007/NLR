import argparse

def bed2dict(file):
    """
    将bed文件转化为以基因名为键的字典
    :param file: bed文件路径
    :return:{"基因名":["染色体",起点,终点]}
    """
    result_dict = {}
    with open(file, 'r') as f:
        line = True
        while line:
            line = f.readline()
            if not line:
                break
            line_list = line.strip().split()
            result_dict[line_list[3]] = line_list[0:3]
    return result_dict


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


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-b","--bed",help='输入bed文件')
    parser.add_argument("-c","--csv",help='输入nlr_summary.csv文件')
    parser.add_argument("-o","--out",help='输出文件路径')
    args = parser.parse_args()

    bed_dict = bed2dict(args.bed)
    # bed文件转化为字典：{"基因名":["染色体",起点,终点]}

    nlr_dict = nlr_summary2dict(args.csv)
    # nlr_summary.csv文件转化为字典：{"基因名":"基因名\t evidence\t structure"}

    # nlrbed_dict记录了nlr基因及其所在染色体、起点和终点：{"基因名":["染色体",起点,终点]}
    nlrbed_dict = {}
    for i in nlr_dict:
        nlrbed_dict[i] = bed_dict[i]

    # 对nlrbed_dict根据染色体，和起点排序，生成nlrbed_items
    # [
    # (基因1,["染色体",起点,终点]),
    # (基因2,["染色体",起点,终点]),...
    # (基因n,["染色体",起点,终点])
    # ]
    nlrbed_items = sorted(nlrbed_dict.items(),key=lambda item:(item[1][0],item[1][1]))

    # 将成簇的序列存储到cluster_dict中：{"基因名":"基因名\t evidence\t structure"}
    # 添加cluster计数功能
    with open(args.out,'w') as f:
        f.write("cluster\tgene\tevidence\tstructure\n")
        cluster_dict = {}
        i = 0
        n = 1
        while i < (len(nlrbed_items)-1):
            if (nlrbed_items[i][1][0] == nlrbed_items[i+1][1][0]) & (int(nlrbed_items[i][1][1])-int(nlrbed_items[i+1][1][2]) < 200000):
                cluster_dict[nlrbed_items[i][0]] = f"{n}\t{nlr_dict[nlrbed_items[i][0]]}"
                cluster_dict[nlrbed_items[i+1][0]] = f"{n}\t{nlr_dict[nlrbed_items[i+1][0]]}"
            elif nlrbed_items[i][1][0] != nlrbed_items[i+1][1][0]:
                n += 1
            i += 1
        for i in cluster_dict:
            f.write(cluster_dict[i])