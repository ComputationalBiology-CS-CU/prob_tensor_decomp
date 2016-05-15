import numpy


if __name__ == "__main__":
    outlier = 142
    targetid = ""
    sid_map = dict()
    with open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_chromosome_22", "r") as f:
        i = 0
        for line in f:
            line = line.strip()
            cols = line.split('\t')
            for j in range(1, len(cols)):
                sid = cols[j][:10]
                if sid not in sid_map:
                    sid_map[sid] = []
                sid_map[sid].append(j-1)

                if j==outlier+1:
                    targetid=sid
            break
            i += 1


    print sid_map[targetid]

    arr = ["" for i in range(1067)]
    with open("../data_processed/sample_subtissue","r+") as f1:
        for lines in f1:
            lines = lines.strip()

            cols = lines.split('\t')
            for i in range(int(cols[1]), int(cols[2])+1):
                arr[i] = cols[0]

    for ele in sid_map[targetid]:
        print arr[ele]
