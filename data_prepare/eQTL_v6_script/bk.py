
    file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_" + str(filter), 'r')
    eQTL_tissue = {}
    while 1:
        line = (file.readline()).strip()
        if not line:
            break

        tissue = (line.split('\t'))[0]
        eQTL_tissue[tissue] = []
    file.close()

    ## sample_list
    file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype", 'r')
    sample_list = (((file.readline()).strip()).split('\t'))[1:]
    file.close()

    ## sample_tissue_rep
    file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_sample_tissue_type", 'r')
    sample_tissue_rep = {}
    while 1:
        line = file.readline()[:-1]
        if not line:
            break

        line = line.split('\t')

        if len(line) < 3:
            print line
            continue

        sample = line[0]
        tissue = line[2]

        sample_tissue_rep[sample] = tissue
    file.close()

    ## fill in the eQTL_tissue rep
    for sample in sample_list:
        tissue = sample_tissue_rep[sample]
        if tissue in eQTL_tissue:
            eQTL_tissue[tissue].append(sample)

    # save the rep
    file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_sample_" + str(filter), 'w')
    for tissue in eQTL_tissue:
        file.write(tissue + '\t')
        for sample in eQTL_tissue[tissue]:
            file.write(sample + '\t')
        file.write('\n')
    file.close()

    ##============ process the rpkm matrix to get eQTL samples ==============
    ## get the sample_rep first
    sample_rep = {}
    file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_sample_" + str(filter), 'r')
    while 1:
        line = (file.readline()).strip()
        if not line:
            break

        line = line.split('\t')[1:]
        for sample in line:
            sample_rep[sample] = 1
    file.close()


    file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype", 'r')
    file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample", 'w')

    # filter all the samples again
    index_rep = {}
    line = (file.readline()).strip()
    line = line.split('\t')
    file1.write(line[0] + '\t')
    for i in range(1, len(line)):
        sample = line[i]
        if sample in sample_rep:
            index_rep[i] = 1
            file1.write(sample + '\t')
    file1.write('\n')
    
    while 1:
        line = (file.readline()).strip()
        if not line:
            break

        line = line.split('\t')
        file1.write(line[0] + '\t')
        for i in range(1, len(line)):
            if i in index_rep:
                file1.write(line[i] + '\t')
        file1.write('\n')

    file.close()
    file1.close()


    ##======================================================================================================
    ##==== remove all the NULL genes as defined (testing for all samples)
    ##==== target: GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null
    ##======================================================================================================
    file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample", 'r')
    file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null", 'w')
    line = file.readline()
    file1.write(line)

    while 1:
        line = (file.readline()).strip()
        if not line:
            break

        # check this gene
        line = line.split('\t')
        if check_null(line[1:]):
            continue
        else:
            file1.write(line[0] + '\t')
            for i in range(1, len(line)):
                file1.write(line[i] + '\t')
            file1.write('\n')

    file.close()
    file1.close()
    
    ##==== remove temp data
    os.system("rm " + "../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample")




    '''
    ##=============================================================================================
    ##==== normalizing all the samples (here we use Log normalize other than the previous Quantile)
    ##==== target: GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize
    ##=============================================================================================
    ## the following is the log normalization

    normal_quantile = 1
    normal_log = 0


    if normal_quantile:
        print "quantile normalization..."

        file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_chromosome_"+str(cid), 'r')
        file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_normalize", 'w')
        line = file.readline()
        file1.write(line)


        gene_list = []
        description_list = []
        rpkm_rep = {}  # hashing gene to L1
        L2 = []

        while 1:
            line = (file.readline()).strip()
            if not line:
                break

            ## sorting; indexing (additive filling L2)
            line = line.split('\t')
            gene = line[0]
            description = line[1]
            gene_list.append(gene)
            description_list.append(description)

            expression = map(lambda x: float(x), line[2:])
            expression = np.array(expression)

            sort = np.argsort(expression)
            ## get the ordering list (or the rank list for this gene for all samples)
            sort1 = []
            for i in range(len(sort)):
                sort1.append(0)
            for i in range(len(sort)):  # the current i is the rank
                index = sort[i]
                sort1[index] = i

            rpkm_rep[gene] = sort1

            if len(L2) == 0:
                for pos in sort:
                    rpkm = expression[pos]
                    L2.append(rpkm)
            else:
                for i in range(len(sort)):
                    pos = sort[i]
                    rpkm = expression[pos]
                    L2[i] += rpkm
        file.close()

        length = len(gene_list)
        for i in range(len(L2)):
            L2[i] = L2[i] * 1.0 / length


        for i in range(len(gene_list)):
            ## two lists:
            ## L1 (value as the re-mapped positions of all original elements; each gene has one such list)
            ## L2 (containing the normalized/averaged value for each index/position)
            gene = gene_list[i]
            description = description_list[i]
            L1 = rpkm_rep[gene]

            file1.write(gene + '\t' + description + '\t')

            for index in L1:
                value = L2[index]
                file1.write(str(value) + '\t')

            file1.write('\n')

        file1.close()

        ##==== remove temp data
        os.system("rm " + "../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null")


    if normal_log:
        print "log normalization..."

        file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null", 'r')
        file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize", 'w')
        line = file.readline()
        file1.write(line)

        while 1:
            line = (file.readline()).strip()
            if not line:
                break

            line = line.split('\t')
            gene = line[0]
            file1.write(gene + '\t')

            rpkm_list = map(lambda x: float(x), line[1:])
            for i in range(len(rpkm_list)):
                rpkm = rpkm_list[i]
                rpkm = math.log(rpkm + 0.1) # NOTE: here is the rule of transformation: shifted logarithm
                file1.write(str(rpkm) + '\t')

            file1.write('\n')

        file.close()
        file1.close()

        ##==== remove temp data
        os.system("rm " + "../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null")



    
    '''



    #### I won't split the samples into their tissues, as I will load the full tensor into memory
    ##======================================================================================================
    ##==== separating all the esamples into their tissues
    ##==== target: expression_by_etissue/tissue_list.txt
    ##==== target: expression_by_etissue/tissue_x.txt
    ##======================================================================================================
    """
    # get the etissue list
    eQTL_tissue = {}
    file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_sample_" + str(filter), 'r')
    while 1:
        line = (file.readline()).strip()
        if not line:
            break

        line = line.split('\t')
        tissue = line[0]
        eQTL_tissue[tissue] = {}
        for i in range(1, len(line)):
            sample = line[i]
            eQTL_tissue[tissue][sample] = 1

    # tissue list
    file = open("../data_processed/expression_by_etissue/tissue_list.txt", 'w')
    tissue_list = []
    count = 0
    for tissue in eQTL_tissue:
        tissue_list.append(tissue)
        count += 1
        file.write(tissue + '\t' + str(count) + '\n')
    file.close()

    # extract esamples into their etissues
    count = 0
    for tissue in tissue_list:
        count += 1
        file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize", 'r')
        file1 = open("../data_processed/expression_by_etissue/tissue_" + str(count) + ".txt", 'w')

        # filter all the samples again
        index_rep = {}
        line = (file.readline()).strip()
        line = line.split('\t')
        file1.write(line[0] + '\t')
        for i in range(1, len(line)):
            sample = line[i]
            if sample in eQTL_tissue[tissue]:
                index_rep[i] = 1
                file1.write(sample + '\t')
        file1.write('\n')
        
        while 1:
            line = (file.readline()).strip()
            if not line:
                break

            line = line.split('\t')
            file1.write(line[0] + '\t')
            for i in range(1, len(line)):
                if i in index_rep:
                    file1.write(line[i] + '\t')
            file1.write('\n')

        file.close()
        file1.close()
    """

