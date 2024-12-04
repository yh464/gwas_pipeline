# stand-alone script to quality control SNPs
f = open('../genqc/exclude_list.txt','w')
chrs = [str(i) for i in range(1,23)] + ['X','XY']
for c in chrs:
    mfi = open(f'/rds/project/rb643/rds-rb643-ukbiobank2/Data_Genetics/ukb_mfi_chr{c}_v3.txt')
    while True:
        line = mfi.readline()
        if len(line) == 0: break
        line = line.replace('\n','').split('\t')
        try:
            line[-1] = float(line[-1])
            line[-3] = float(line[-3])
            if line[-1] < 0.4 or line[-3] < 0.001: print(line[1], file = f)
        except:
            print(line[1],file = f)