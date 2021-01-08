
print("hello")

import re #Import regular expressions parser

#Import .agp file, map, and designate output file
agp_file_path = input("Specify path to agp file: ")
map_file_path = input("Specify path to map file: ")
outfile_path = input("Specify path to output file: ")

#AGP
#E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/Lep-Map3/Final_Map_F2_F3_F5/AGPs_final/full_map_updated_correct.agp
#E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/Lep-Map3/Final_Map_F2_F3_F5/Final_map_joined_LGs/Final_AGPs_modified/full_genome.agp #Final AGP with joined LGs 6-8-2020

#Map file
#E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/QTL2/Final_map/ipomoea_gmap.csv #F2 map
#E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/Lep-Map3/Final_Map_F2_F3_F5/Final_map_joined_LGs/Converted_linkage_maps/F3_map.csv  #F3map
#E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/Lep-Map3/Final_Map_F2_F3_F5/Final_map_joined_LGs/Converted_linkage_maps/F5_pre_conversion.csv  #F3map

#Out file
#E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/Lep-Map3/Final_Map_F2_F3_F5/Final_map_joined_LGs/Converted_linkage_maps/converted_final_F2_map_6_8_2020.csv
#E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/Lep-Map3/Final_Map_F2_F3_F5/Final_map_joined_LGs/Converted_linkage_maps/converted_final_F3_map_6_8_2020.csv
#E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/Lep-Map3/Final_Map_F2_F3_F5/Final_map_joined_LGs/Converted_linkage_maps/converted_final_F5_map_6_8_2020.csv


#This block of code reads in the .agp file as a dictionary where each scaffold is a key that connects to position ranges in the scaffold and in the chromosome
agp={}

#print(marker_file)
linecounter = 0
with open(agp_file_path, "r") as f:
    for line in f:
        if line[0] != "#":
            print(line)
            line = line.strip()
            values = re.split('\t', line)
            print(values)
            #print(values[5])
            #print(values[5][0])
            if values[5][0] == "S":
                chrrLG = values[0]
                print(chrrLG)
                LG_pos = (int(values[1]), int(values[2]))
                print(LG_pos)
                scaff_pos = (int(values[6]), int(values[7]))
                print(scaff_pos)
                orientation=values[8]
                print(orientation)
                scaffold=int(re.split('\_', values[5])[1])
                print(scaffold)
                #agp[scaffold].append(scaff_pos, LG_pos, orientation)
                #Stuck on making list of tuples that is not class tuple
                if scaffold not in agp.keys():
                    #agp[scaffold]={}
                    #agp[scaffold]=list(scaff_pos, LG_pos, orientation)
                    agp[scaffold]=[]
                    print("key")
                    print(agp)
                    print(agp[scaffold])
                    agp[scaffold].append((scaff_pos, LG_pos, orientation, chrrLG))
                    #agp[scaffold][scaff_pos]=(LG_pos, orientation)
                else:
                    agp[scaffold].append((scaff_pos, LG_pos, orientation, chrrLG))
        linecounter += 1
        print("line count ", linecounter)

print(agp)
print(agp.keys())
values=[]

#This block of code reads the map file and for each scaffold-position combination in the map file, finds the corresponding scaffold and position in the .agp file and converts the position on the scaffold to the position on the chromosome, taking scaffold orientation into account
linecounter = 0
outlines=[]
with open(map_file_path, "r") as f:
    for line in f:
        if linecounter>0:
        #if line[0] != "#":
           # print(line)
            line = line.strip()
            values = re.split(',', line)
            print(values)
            print(values[0])
            converted_values = []
            to_convert = values[0]
            print(to_convert)
        #print(linecounter)
                if to_convert[0]=="S":
                    scaffold_to_convert = int(re.split('\_', to_convert)[1])
                    position_to_convert = int(re.split('\_', to_convert)[2])
                    print("scaffold", scaffold_to_convert)
                    print(position_to_convert)
                    print(agp[scaffold_to_convert])
        #linecounter+=1
                    if scaffold_to_convert in agp.keys():
                        print("FOUND")
                        for i in agp[scaffold_to_convert]:
                           print("scaffold entry ", i)
    #        linecounter += 1
                        if (i[0][0] < position_to_convert < i[0][1]) == True:
                                print("YIPPEE")
                                print(position_to_convert)
                                print(i[2])
                                LG = i[3]
                                if (i[2] == '+'):
                                    print("Forward")
                                    new_position = (((position_to_convert - i[0][0])) + (i[1][0]))
                                elif (i[2] == '-'):
                                    print("reverse")
                                    new_position = ((i[0][1] - position_to_convert) + i[1][0])
                                print(new_position)
                else:
                    new_position = "NA"
                converted_values.append(values[0])
                converted_values.append(values[1])
                converted_values.append(values[2])
                converted_values.append(i[3])
                converted_values.append(new_position)
                outlines.append(converted_values)
            linecounter += 1

print(outlines)



#This block of code creates the output file

#f = open("AGP_converted.csv", "w+")
f = open(outfile_path, "w+")

#f.write("Marker_number, Scaffold, Position_on_scaffold, LG_in_ASMAP, CM_position, Position_on_LG, LG,\n")
f.write("Marker, LG, CM, AGP_LG, Position_on_LG\n")
#6830,17266,33348559,5,1.818466791,5650015,L.1,
#print(outfile[0])
for i in range(0,len(outlines)):
#    print((outfile[i]))
#    print(str((outfile[i])))
    for j in (outlines[i]):
        f.write(str(j))
        f.write(",")
        print(str(j))
    f.write("\n")
f.close()
