#!user/bin/python

#SRR1916153
#SRR1916156
#SRR1916155
#SRR1916154
#SRR1916152
# fastq-dump accession
# fastqc accession.fastq

import sys
acc_file = sys.argv[1]
with open(acc_file) as acc_file_handle:
    for line in acc_file_handle:
        #print(line)
        line = line.strip()
        cmd_str = "fastq-dump {srr_acc}".format(srr_acc = line)
        print(cmd_str)
