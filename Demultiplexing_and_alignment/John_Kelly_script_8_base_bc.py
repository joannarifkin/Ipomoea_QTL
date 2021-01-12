'''
Created on Nov 28, 2015

@author: Joanna
'''

if __name__ == '__main__':
    pass

#!/usr/bin/env python

import sys
import re

# parser for paired end data with diagnostics 
barcode_file= open("complete_barcodes_file.txt","rU") #"key2014r.txt"

fq2 = sys.argv[2]#"testR2_section2_copied.txt"#"testR2_copied.txt"#"testR2_section2_copied.txt"#"testR2.fastq"#sys.argv[2]
fq1 = sys.argv[1]#"testR1_section2_copied.txt"#"testR1_copied.txt"#"testR1_section2_copied.txt"#testR1.fastq"#sys.argv[1]
index=sys.argv[3]#"ATCACG"#"CGATGT"#"ATCACG"#sys.argv[3]

sumbc=open(index+"bcx.txt","w")

barcodes_dict = dict()
mainbc=dict()
for idx,row in enumerate(barcode_file):
    line = row.replace('\n', '').replace('\r', '').split('\t')
    if line[0]==index:
        barcodes_dict[line[1]] = {'name':line[2], 'file_name':line[2]+".fq"}
        mainbc[line[1]]=1

allbc=dict()
bases=["A","C","G","T"]
for k1 in range(4):
    for k2 in range(4):
        for k3 in range(4):
            for k4 in range(4):
                for k5 in range(4):
                    for k6 in range(4):
                        for k7 in range(4):
                            for k8 in range(4):
                                key=bases[k1]+bases[k2]+bases[k3]+bases[k4]+bases[k5]+bases[k6]+bases[k7]+bases[k8]
                                allbc[key]=[0,0]
                    

r1_file = open(fq1, 'r')
r2_file = open(fq2, 'r')
current_file1 = ''
current_file2 = ''

matches=[0,0,0]
for idx, row in enumerate(r1_file):
    line = row.replace('\n', '').replace('\r', '')
    row2 = r2_file.readline()
    line2= row2.replace('\n', '').replace('\r', '')

    if idx%4 == 0:
        name_line = row
        name_line2 = row2                
    elif idx%4 == 1:
        catok=[0,0]
        try:
            zipa=mainbc[row[:8]]
            catok[0]=1
        except KeyError:
            pass
        try:
            zipa=mainbc[row2[:8]]
            catok[1]=1
        except KeyError:
            pass

        if catok[0]==1 and catok[1]==1 and row[:8] == row2[:8]:
            matches[0]+=1
            current_file1 = open("r1."+barcodes_dict[row[:8]]['file_name'], 'a')
            current_file1.write(name_line)
            current_file1.write(row[8:])
            current_file2 = open("r2."+barcodes_dict[row[:8]]['file_name'], 'a')
            current_file2.write(name_line2)
            current_file2.write(row2[8:])                        
        elif (catok[0]==0 and catok[1]==0) or (catok[0]==1 and catok[1]==1 and row[:8] != row2[:8]):
            matches[2]+=1
            current_file1 = open(index+'misc.barcodes1.fq', 'a')
            current_file1.write(name_line)
            current_file1.write(row[:])       #Altering to keep the entire read rather than just the part after the barcode                 
            current_file2 = open(index+'misc.barcodes2.fq', 'a')
            current_file2.write(name_line2)
            current_file2.write(row2[:]) 
        else:
            matches[1]+=1
            if catok[0]==1:
                current_file1 = open("r1."+barcodes_dict[row[:8]]['file_name'], 'a')
                current_file1.write(name_line)
                current_file1.write(row[8:])
                current_file2 = open("r2."+barcodes_dict[row[:8]]['file_name'], 'a')
                current_file2.write(name_line2)
                current_file2.write(row2[8:])
            else:
                current_file1 = open("r1."+barcodes_dict[row2[:8]]['file_name'], 'a')
                current_file1.write(name_line)
                current_file1.write(row[8:])
                current_file2 = open("r2."+barcodes_dict[row2[:8]]['file_name'], 'a')
                current_file2.write(name_line2)
                current_file2.write(row2[8:])                                
        try:
            allbc[row[:8]][0]+=1
            allbc[row2[:8]][1]+=1
        except KeyError:
            pass


    elif idx%4 == 2:
        current_file1.write(row)
        current_file2.write(row2)                
    else:
        current_file1.write(row[8:])
        current_file2.write(row2[8:])


z=allbc.keys()
for zz in z:
    if allbc[zz][0]>0 or allbc[zz][1]>0:
        try:
            xt=mainbc[zz]
            sumbc.write(zz+'\t'+str(allbc[zz][0])+'\t'+str(allbc[zz][1])+'\t1\n')
        except:
            sumbc.write(zz+'\t'+str(allbc[zz][0])+'\t'+str(allbc[zz][1])+'\t0\n') 

sumbc.write('000\t'+str(matches[0])+'\t'+str(matches[1])+'\t'+str(matches[2])+'\n')


