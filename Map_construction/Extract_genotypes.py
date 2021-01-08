
print("hello")

import re #Import regular expressions parser

# make smaller post call file in command line
# cut -f 1-4 post_call_sample_with_header.txt > post_call_sample_with_header_LPR_CAA_only.txt

#genotype likelihoods (posteriors) for each 10 possible SNP genotypes AA, AC, AG, AT, CC, CG, CT, GG, GT and TT.

post_call_file_path = input("Specify path to post.call file: ")
#post_call_sample_with_header_LPR_CAA_only.txt
#mini_post_call_sample_with_header.txt
#post_call_sample_with_header.txt
#test_NA_genotypes.txt
out_file_path = input("Specify path to output file")
#converted_post_call_sample_with_header_LPR_CAA_only.txt
#converted_post_call_sample.txt
#mini_converted_post_call_sample.txt

linenumber=0
header=[]
outfile=[]
with open(post_call_file_path, "r") as f:
    for line in f:
        if linenumber<7:
            header.append(line)
        if linenumber>=7:
            print(line)
            values=[]
            values = re.split('\t', line)
            print(values)
            scaffold=values[0]
            print(scaffold)
            print(len(values))
            position=values[1]
            print(position)
            in_genotypes=values[2:len(values)]
            print(in_genotypes)
            out_genotypes=[]
            out_genotypes.append(scaffold)
            out_genotypes.append(position)
            for gen in in_genotypes:
                #if gen == "1 1 1 1 1 1 1 1 1 1":
                 #   out_genotypes.append("NA")
                #elif gen == "1.0":
                  #  out_genotypes.append("NA")
                if gen == "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0":
                #elif gen == "1.0\s1.0\s1.0\s1.0\s1.0\s1.0\s1.0\s1.0":
                    out_genotypes.append("NA")
                    print("NA")
                elif gen == "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n":
                    out_genotypes.append("NA")
                    print("NA")
                else:
                    genvalues=re.split('\s', gen)
                    print(genvalues)
                    print("genotypes n")
                    print(len(genvalues))
                    if genvalues[0]=="1.0":
                        out_genotypes.append("AA")
                    if genvalues[1] == "1.0":
                        out_genotypes.append("AC")
                    if genvalues[2] == "1.0":
                        out_genotypes.append("AG")
                    if genvalues[3] == "1.0":
                        out_genotypes.append("AT")
                    if genvalues[4] == "1.0":
                        out_genotypes.append("CC")
                    if genvalues[5] == "1.0":
                        out_genotypes.append("CG")
                    if genvalues[6] == "1.0":
                        out_genotypes.append("CT")
                    if genvalues[7] == "1.0":
                        out_genotypes.append("GG")
                    if genvalues[8] == "1.0":
                        out_genotypes.append("GT")
                    if genvalues[9] == "1.0":
                        out_genotypes.append("TT")
                print(gen)
                #print(out_genotypes)
            outfile.append(out_genotypes)
        linenumber+=1
        print(linenumber)
    outfile.append(header)

in_genotypes[-1]

print(outfile)
print(values[3])
print(header)
print(header[2])
print(str.split(header[2],"\t"))
print(len(str.split(header[2],"\t")))
print(outfile[0])

str.split(header[2],"\t")[-1]


f = open(out_file_path, "w+")
#f.write("Marker_number, Scaffold, Position_on_scaffold, LG_in_ASMAP, CM_position, Position_on_LG, LG,\n")
#6830,17266,33348559,5,1.818466791,5650015,L.1,
#print(outfile[0])
#for i in range(0, len(header)):
for j in range(len(str.split(header[2],"\t"))-1):
    print(j)
#    f.write(str(j))
    f.write(str(str.split(header[2],"\t")[j]))
    f.write(",")
f.write(str(str.split(header[2],"\t")[-1]))
#f.close()
#f.write("\n")
for i in range(0,len(outfile)):
#    print((outfile[i]))
#    print(str((outfile[i])))
    for j in range(len(outfile[i])-1):
        f.write(str(outfile[i][j]))
        f.write(",")
        print(str(j))
    f.write(str(outfile[i][-1]))
    f.write("\n")
f.close()