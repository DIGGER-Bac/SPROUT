# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 09:20:18 2023

@author: Benedikt Wiebach
Email: bndkt.wbch@gmail.com
"""

##############################################################################

# imports
import subprocess
import os
import sys
import statistics as s
import time
import shutil

path = os.getcwd()                                                             # change working directory
try:
    os.mkdir(path+"/Data/processing")
except OSError:
    print("")
#copy perl-script in data folder
shutil.copy("/usr/local/share/ViennaRNA/bin/relplot.pl",path+"/Data/processing/relplot.pl") # copy script of ViennaRNA for positional entropy tothe right position
os.chdir(path +'/Data')                                                        # change working directory



filenames = open("files.txt","w")                                              # read out filenames from program folder
subprocess.run(["ls"], stdout=filenames, shell = True)          
filenames.close()

filenames = open("files.txt","r")
files = filenames.readlines()
filenames.close()

flag_stop = False
while (flag_stop == False):
    for file in files:                                                         # only ".fasta" datatype is interesting
        if ".fasta" in file:
            flag_stop = True
            continue
        else:
            files.remove(file)
            flag_stop = False
            if files ==[]:                                                     # stop programm if there are no files to be analyzed
                flag_stop = True
                print("There are no sequences in .fasta format to be analyzed in the folder "+path+"\Input")
                sys.exit()
            break
i = -1
for file in files:
    i = i + 1
    f = open(file.strip(),"r")
    temp = f.readlines()
    z = -1
    for line in temp:                                                          # cut out only the sequence
        z = z + 1
        if ">" in line:
            path = os.getcwd()                                                             # change working directory
            os.chdir(path +'/processing')                                                        # change working directory
            new_filename = open(line.split(">")[1].strip()+".fasta","w") 
            new_filename.write(line)
            new_filename.write(temp[z+1])
            new_filename.close()
            os.chdir('..')
        else:
            continue
    f.close()




os.remove("files.txt")
path = os.getcwd()                                                             # change working directory
os.chdir(path +'/processing')                                                        # change working directory

#"""
##############################################################################
# preparing the sequences
##############################################################################

filenames = open("files.txt","w")                                              # read out filenames from program folder
subprocess.run(["ls"], stdout=filenames, shell = True)          
filenames.close()

filenames = open("files.txt","r")
files = filenames.readlines()
filenames.close()

flag_stop = False
while (flag_stop == False):
    for file in files:                                                         # only ".fasta" datatype is interesting
        if ".fasta" in file:
            flag_stop = True
            continue
        else:
            files.remove(file)
            flag_stop = False
            if files ==[]:                                                     # stop programm if there are no files to be analyzed
                flag_stop = True
                print("There are no sequences in .fasta format to be analyzed in the folder "+path+"\Data")
                sys.exit()
            break

Seq = [0 for i in range(len(files))]                                           # create an array with the size of the amount of sequences
i = 0
for file in files:
    f = open(file.strip(),"r")
    temp = f.readlines()
    for line in temp:                                                          # cut out only the sequence
        if ">" in temp:
            continue
        else:
            Seq[i] = line
    f.close()
    i = i+1

# prediction of the stem loop and its probability
# gathering of needed data
 
Seq_Prob_names = [0 for i in range(len(files))]  
Seq_2DStruct = []
i = 0
for file in files:    
    out_RNAfold = open("RNAfold.txt","w")   
    subprocess.run(["RNAfold", "-p", "--MEA", file.strip()], stdout=out_RNAfold, shell=False ) # execute RNAfold from the ViennaRNA package         
    out_RNAfold.close()
    print("RNAfold Seq ", i+1)
    out_RNAfold = open("RNAfold.txt","r")                                      #save prediction for secondary structure
    temp2 = out_RNAfold.readlines()
    Seq_2DStruct.append(temp2[2])
    out_RNAfold.close()

    tmp_file = file.strip().split(".")                                         # produce data names needed for relplot.pl
    tmp_str1 = tmp_file[0] + "_ss.ps"
    tmp_str2 = tmp_file[0] + "_dp.ps"
    tmp_str3 = tmp_file[0] + "_rss.ps"
    handle_rss = open(tmp_file[0] + "_rss.ps","w")
    
    Seq_Prob_names[i] = tmp_str3
    i = i+1
    
    try:
        #print(tmp_file[0] + "_ss.ps")
        handle_ss = open(tmp_file[0] + "_ss.ps","r")
        handle_ss.close()
    except FileNotFoundError:                                                  # stop program if there is an error with the file naming
        print("The file " + file.strip() + " is named wrong.")
        print(" The name of the file is not identical with the name of the sequence.")
        print(" Consider the naming rules for input files in the ReadMe and try again.")
        try:
            os.remove(tmp_str1)
        except FileNotFoundError:                                                        # delete temporary variables
            t = 0
        try:
            os.remove(tmp_str2)
        except FileNotFoundError:
            t = 0
        try:
            os.remove(tmp_str3)
        except FileNotFoundError:
            t = 0
        sys.exit()  
    
    subprocess.run([ "perl", "relplot.pl", tmp_str1, tmp_str2], stdout=handle_rss, shell=False ) # use script for positional entropy
    handle_rss.close()
    
    os.remove(tmp_str1)                                                        # delete temporary variables
    os.remove(tmp_str2)

os.remove("relplot.pl")

Seq_Prob = [[] for i in range(len(files))]                                     # datastructure is a 2D array
i = -1
for file in Seq_Prob_names:                                                    # same order of sequences as for structure prediction
    i = i+1
    handle = open(file,"r")
    lines = handle.readlines()
    flag_prob = False   
    for line in lines:
        if "/S [" in line:                                                     # start condition for reading in data
            flag_prob = True
            continue
        else:
            if flag_prob == True:
                    if "] def" in line:                                        # stop condition for reading in data
                        flag_prob = False
                        continue
                    else:
                        tmp3 = line
                        Seq_Prob[i].append(float(tmp3.strip()))                               # add to list
                        continue
            else:    
                continue
    handle.close()
time.sleep(1)
for file in Seq_Prob_names:                                                    # delete no longer needed data
    os.remove(file)
for file in files:                                                    # delete no longer needed data
    os.remove(file.strip())    
os.remove("files.txt")
os.remove("RNAfold.txt")

i = 0                   
for sequence in Seq:                                                           # invert the sequence
    Seq[i] = sequence[::-1]
    i=i+1
    
i = 0
for sequence in Seq_2DStruct:                                                   # invert secondary structure sequence
    tmp4 = sequence.split(" ")
    Seq_2DStruct[i] = tmp4[0][::-1]                                             # WARNING: brackets are now shown the opposite way!
    i=i+1
    
i=0  
for sequence in Seq_Prob:                                                      # invert structure entropy
    Seq_Prob[i] = list(reversed(Seq_Prob[i]))              
    i=i+1    

##############################################################################
# collect data from each sequence
##############################################################################
toDelete = []                                                                  # array containing the indexes of unusable datasets

# read in variables from config file
##############################################################################
# (stop program if there is an error with one of the variables)
os.chdir('..')                                                                 # change working directory
os.chdir('..')

configuration_file = open("config.yml","r")
configuration = configuration_file.readlines()
configuration_file.close

for line in configuration:
    
    if "#" in line:
        continue
    
    # sequence reading
    if "Skip_noU" in line:
        var = line.split("=")
        if var[1].strip().isdigit():
            Skip_noU = int(var[1].strip())
        else:
            print("'Skip_noU' is not defined as an integer in the configuration file. Please change it and start the program again.")
            sys.exit()
        
    if "threshold_Prob" in line:
        var = line.split("=")
        try:
            threshold_Prob = float(var[1].strip())
        except (ValueError,IndexError):
            print("'threshold_Prob' is not defined as a float in the configuration file. Please change it and start the program again.")
            sys.exit()

    if "threshold_StemLoopError" in line:
        var = line.split("=")
        if var[1].strip().isdigit():
            threshold_StemLoopError = int(var[1].strip())     
        else:
            print("'threshold_StemLoopError' is not defined as an integer in the configuration file. Please change it and start the program again.")
            sys.exit()
            
    if "max_gap_to_UA" in line:
        var = line.split("=")
        if var[1].strip().isdigit():
            max_gap_to_UA = int(var[1].strip())
            max_gap_to_UA = max_gap_to_UA + 1        
        else:
            print("'max_gap_to_UA' is not defined as an integer in the configuration file. Please change it and start the program again.")
            sys.exit()
            
    if "window_size_UA_min" in line:
        var = line.split("=")
        if var[1].strip().isdigit():
            window_size_UA_min = int(var[1].strip())        
        else:
            print("'window_size_UA_min' is not defined as an integer in the configuration file. Please change it and start the program again.")
            sys.exit()
            
    if "window_size_UA_max" in line:
        var = line.split("=")
        if var[1].strip().isdigit():
            window_size_UA_max = int(var[1].strip())      
        else:
            print("'window_size_UA_max' is not defined as an integer in the configuration file. Please change it and start the program again.")
            sys.exit()
            
    # sequence generating    
               
    if "amount_sequences" in line:
        var = line.split("=")
        if var[1].strip().isdigit():
            amount_sequences = int(var[1].strip())        
        else:
            print("'amount_sequences' is not defined as an integer in the configuration file. Please change it and start the program again.")
            sys.exit()
             

#"""

flag_replace_sRNA_generate_config = True

path = os.getcwd()                                                             # change working directory
os.chdir(path +'/Data')

# create datastructure
##############################################################################
# for each sequence:   | 0 start-pos | 1 end-pos | 2 length | 3 amount U or  loop-size | 4 amount A
                                                                      # -> loop-size is length of the loop
data_polyU = [[0,0,0,0] for i in range(len(files))]                            # create datastructure for analysis
data_stemL = [[0,0,0,0] for i in range(len(files))]
data_richUA = [[0,0,0,0,0] for i in range(len(files))]

# recognization Poly-U
##############################################################################

i = 0                                                                          # sequence number
for sequence in Seq:                                                           # Poly-U
    flag_start_U = False
    flag_NoU = False
    NoU_count = 0
    j=1

    for base in sequence:
        if NoU_count == Skip_noU and flag_NoU == False:
            flag_NoU = True                                                    # NoU                                                                              
        if base == 'U' and flag_start_U == False:
            data_polyU[i][0] = j
            data_polyU[i][3] = data_polyU[i][3] + 1
            flag_start_U = True
            flag_NoU = False
            NoU_count = 0                                                      # start                                                                              
        elif base == 'U' and flag_start_U == True:
            data_polyU[i][3] = data_polyU[i][3] + 1
            data_polyU[i][1] = j
            flag_NoU = False
            NoU_count = 0                                                      # U                                                                               
        if base != 'U' and flag_NoU == False:
            NoU_count = NoU_count + 1                                          # no U                                                                               
            j = j + 1
            continue
        elif base != 'U' and flag_NoU == True:
            if flag_start_U == False:
                data_polyU[i][2] = 0
                toDelete.append(i)
            if flag_start_U == True:
                data_polyU[i][2] = data_polyU[i][1] - data_polyU[i][0] + 1     # stop                                                                       
            break
        j = j + 1
    i = i+1

# recognization stemloop
##############################################################################

i = 0                                                                          # sequence number
flag_noStemLoopFound = False
for sequence in Seq_2DStruct:                                                  # stem loop
    flag_start_Loop = False
    flag_in_Loop = False
    counter_stemL_error = 0
    count_brackets_l = 0
    count_false_brackets_l = 0
    count_brackets_r = 0
    count_false_brackets_r = 0
    j=1

    for base in sequence:                                                      # WARNING: opening and closing brackets are inverted
        if base == "." and flag_in_Loop == False:
            j = j + 1
            continue
        
        elif base == "." and flag_in_Loop == True and count_brackets_l-count_brackets_r > 2:
            data_stemL[i][3] = data_stemL[i][3] + 1
            j = j + 1
            counter_stemL_error = 0
            continue
        elif base == "." and flag_in_Loop == True:
            counter_stemL_error = counter_stemL_error + 1
            
            if count_brackets_l + count_false_brackets_l - count_brackets_r - count_false_brackets_r < 2: # intercept possibility that opening and closing brackets have different probabilities
                data_stemL[i][1] = j - counter_stemL_error
                data_stemL[i][2] = data_stemL[i][1] - data_stemL[i][0] + 1
                break
            if counter_stemL_error > threshold_StemLoopError:
                flag_noStemLoopFound = True                                    # delete dataset, if no stem loop is found
                data_stemL[i][1] = 0
                data_stemL[i][1] = 0
                data_stemL[i][2] = 0
                data_stemL[i][3] = 0
                toDelete.append(i)
                break
            else:
                j = j + 1
                continue
        
        if base == ")" and flag_in_Loop == False:                              # start stemloop
            if Seq_Prob[i][j-1] <= threshold_Prob:
                count_brackets_l = count_brackets_l + 1
                data_stemL[i][0] = j
                flag_in_Loop = True
            else:
                count_false_brackets_l = count_false_brackets_l +1
        
        elif base == ")" and flag_in_Loop == True:                             # extend stemloop
            if Seq_Prob[i][j-1] <= threshold_Prob:
                count_brackets_l = count_brackets_l + 1
                counter_stemL_error = 0
            else:
                count_false_brackets_l = count_false_brackets_l +1

        elif base == "(":                                                      # close stemloop
            if Seq_Prob[i][j-1] <= threshold_Prob:
                count_brackets_r = count_brackets_r + 1
                counter_stemL_error = 0
            else:
                count_false_brackets_r = count_false_brackets_r +1
       
        if count_brackets_r == count_brackets_l and flag_in_Loop == True:      # stop if the amount of opening and closing brackets is the same
            data_stemL[i][1] = j
            data_stemL[i][2] = data_stemL[i][1] - data_stemL[i][0] + 1
            data_stemL[i][3] = data_stemL[i][3] + count_false_brackets_l + count_false_brackets_r
            break
        if count_brackets_l > 4 and count_brackets_l - count_brackets_r < 3 and Seq_2DStruct[i][j+1] == ")" and Seq_Prob[i][j-1] > threshold_Prob: # intercept possibility that there is a second stemloop direclty after closing of the first one
            data_stemL[i][1] = j - counter_stemL_error                                                       # include error on the end of the stemloop
            data_stemL[i][2] = data_stemL[i][1] - data_stemL[i][0] + 1                                       # if there is an overlap: stemloop ends earlier
            data_stemL[i][3] = data_stemL[i][3] + count_false_brackets_l + count_false_brackets_r
            break
        j=j+1
    i=i+1        


# test if poly-U and stemloop overlap
##############################################################################

i = -1
flag_overlap = False
for sequence in Seq_2DStruct:
    i = i + 1
    
    if data_stemL[i][0] == 0:                                                  # skip if no stemloop was found    
        continue
    if data_polyU[i][1] == 0:                                                  # skif if no poly-U was found    
        continue
    if data_stemL[i][0] > data_polyU[i][1]:                                    # no overlap
        continue
    
    elif data_stemL[i][0] <= data_polyU[i][1]:     
        flag_overlap = True
        while(flag_overlap == True):
            if data_stemL[i][0] == data_polyU[i][1]:                           # if size of overlap is one
                if Seq_Prob[i][data_stemL[i][0]-1] < (threshold_Prob/2):       # if stemloop has very high probability
                    data_polyU[i][1] = data_polyU[i][1] - 1                    # shorten poly-U
                    data_polyU[i][2] = data_polyU[i][2] - 1
                    if Seq[i][data_stemL[i][0]-1] == "U":                      # only decrease amount of U's, if there is a U on that position
                        data_polyU[i][3] = data_polyU[i][3] - 1
                    flag_overlap = False
                    continue
                else:                                                          # otherwise shorten stemloop
                    data_stemL[i][0] = data_stemL[i][0] + 1 
                    data_stemL[i][1] = data_stemL[i][1] - 1
                    data_stemL[i][2] = data_stemL[i][2] - 2 
                    flag_overlap = False
                    continue
            if data_stemL[i][0] < data_polyU[i][1]:                            # if size of overlap is bigger than one
                if Seq_Prob[i][data_stemL[i][0]-1] < (threshold_Prob/2):       # if stemloop has very high probability
                    data_polyU[i][1] = data_polyU[i][1] - 1                    # shorten poly-U
                    data_polyU[i][2] = data_polyU[i][2] - 1
                    #print(Seq[i][data_stemL[i][0]-1], " position: ",data_polyU[i][0]-1+1)
                    if Seq[i][data_polyU[i][1]+1-1] == "U":                      # only decrease amount of U's, if there is a U on that position
                        data_polyU[i][3] = data_polyU[i][3] - 1
                    continue
                else:                                                          # otherwise shorten stemloop
                    data_stemL[i][0] = data_stemL[i][0] + 1 
                    data_stemL[i][1] = data_stemL[i][1] - 1
                    data_stemL[i][2] = data_stemL[i][2] - 2 
                    continue

# recognization UA-rich sequence
##############################################################################

i = -1                                                                         # UA-rich sequence
for sequence in Seq:
    i = i + 1
    window = window_size_UA_min + 1
    flag_windowsize_reached = False
    tempUA = 0
    if data_stemL[i][1] == 0:                                                  # skip if no stemloop was found    
        continue

    while(flag_windowsize_reached == False):
        flag_tempUA_max = False
        if window > window_size_UA_max:
            flag_windowsize_reached = True
            continue
        for gap in range(max_gap_to_UA):
            counter_A = 0
            counter_U = 0
            k = -1 + data_stemL[i][1] + gap                                    # start window after stemloop has ended
            for base in range(window):
                k = k + 1
                if Seq[i][k] == "U":
                    counter_U = counter_U + 1
                elif Seq[i][k] == "A":
                    counter_A = counter_A + 1
            if (counter_U+counter_A)/window >= tempUA:
                if flag_tempUA_max == False:
                    data_richUA[i][0] = data_stemL[i][1] + gap + 1                 # start
                    data_richUA[i][1] = data_stemL[i][1] + gap + window            # end
                    data_richUA[i][2] = data_richUA[i][1] - data_richUA[i][0] + 1  # length
                    data_richUA[i][3] = counter_U                                  # U
                    data_richUA[i][4] = counter_A                                  # A
                    tempUA = (counter_U+counter_A)/window 
                    if tempUA == 1:
                        flag_tempUA_max = True
        window = window + 1

    if tempUA < 0.5:                                                           # if no sufficiently good sequence was found: decrease window size to minimal size
        window = window_size_UA_min
        for gap in range(max_gap_to_UA):
            counter_A = 0
            counter_U = 0
            k = -1 + data_stemL[i][1] + gap                                    # start window after stemloop has ended
            for base in range(window):
                k = k + 1
                if Seq[i][k] == "U":
                    counter_U = counter_U + 1
                elif Seq[i][k] == "A":
                    counter_A = counter_A + 1
            if (counter_U+counter_A)/window >= tempUA:
                data_richUA[i][0] = data_stemL[i][1] + gap + 1                 # start
                data_richUA[i][1] = data_stemL[i][1] + gap + window            # end
                data_richUA[i][2] = data_richUA[i][1] - data_richUA[i][0] + 1  # length
                data_richUA[i][3] = counter_U                                  # U
                data_richUA[i][4] = counter_A                                  # A
                tempUA = (counter_U+counter_A)/window
    if tempUA == 0:
        toDelete.append(i)
    
# delete sequences where not all motives were found
##############################################################################

# toDelete.append(i)  <- command to exclude sequence index

for elements in toDelete:                                                      # avoid double entries
    while(toDelete.count(elements)>1):
        toDelete.remove(elements)

toDelete.sort(reverse = True)                                                  # go through list descending to not change indizes in between

for number in toDelete:                                                        # delete sequences where not all motives were found
    print("sequence ", files [number], " has not all motives. It will not be added to the consensus.")
    if Seq[number] == "":
        print("There is no sequence!")
    del Seq[number]
    del Seq_2DStruct[number]
    del Seq_Prob[number]
    del Seq_Prob_names[number]
    del files[number]
    del data_polyU[number]
    del data_stemL[number]
    del data_richUA[number]

# poly-U
##############################################################################

med_polyU = []
U_percentage_polyU = 0
median_polyU = 0
min_polyU = 0
max_polyU = 0

for entry in data_polyU:
    med_polyU.append(entry[2])
    U_percentage_polyU =  U_percentage_polyU + entry[3]/entry[2]
U_percentage_polyU = U_percentage_polyU / len(data_polyU)
median_polyU = s.median(med_polyU)
min_polyU = min(med_polyU)
max_polyU = max(med_polyU)

print("########################################")
print("poly-U motif:")
print("--------------------")
print("median length: ",median_polyU)
print("average share U: ",U_percentage_polyU)
print("minimal length: ",min_polyU)
print("maximal length: ",max_polyU)
print("########################################")

# stemloop
##############################################################################

med_stemL = []
med_loop = 0

median_stemL = 0
min_stemL = 0
max_stemL = 0
free_bp = 100

for entry in data_stemL:
    med_stemL.append(entry[2])
    med_loop = med_loop + entry[3]

free_bp = med_loop / len(data_stemL)
median_stemL = s.median(med_stemL)
min_stemL = min(med_stemL)
max_stemL = max(med_stemL)

print("stemloop motif:")
print("--------------------")
print("minimal length: ", min_stemL)
print("average length loop  (amount unpaired bases): ",free_bp)
print("average length stem (amount basepairs): ",int(((min_stemL - free_bp)/2)))
print("median length: ",median_stemL)
print("maximal length: ",max_stemL)
print("########################################")

# UA-rich sequence
##############################################################################

med_richUA = []
med_gap =[]
U_percentage = 0
A_percentage = 0
median_richUA = 0
min_richUA = 0
max_richUA = 0
min_gap = 0
max_gap = 0
i = -1
for entry in data_richUA:
    i = i + 1
    med_richUA.append(entry[2])
    med_gap.append(entry[0]-data_stemL[i][1]-1)
    U_percentage =  U_percentage + entry[3]/entry[2]
    A_percentage =  A_percentage + entry[4]/entry[2]
U_percentage = U_percentage / len(data_richUA)
A_percentage = A_percentage / len(data_richUA)

median_richUA = s.median(med_richUA)
min_richUA = min(med_richUA)
max_richUA = max(med_richUA)

median_gap = s.median(med_gap)
min_gap = min(med_gap)
max_gap = max(med_gap)

print("UA-rich motif:")
print("--------------------")
print("median length: ",median_richUA)
print("average share U: ",U_percentage)
print("average share A: ",A_percentage)
print("minimal length: ",min_richUA)
print("maximal length: ",max_richUA)
print("--------------------")
print("gap size to stemloop:")
print("--------------------")
print("median:       ",median_gap)
print("minimal:        ",min_gap)
print("maximal:        ",max_gap)
print("########################################")

##############################################################################
# overwrite config file for sRNA_generate 
##############################################################################

if flag_replace_sRNA_generate_config == True:
    path = os.getcwd() 
    shutil.rmtree(path+"/processing")
    os.chdir('..')
    path = os.getcwd()
    os.chdir(path +"/Output")
    
    new_config = open("consensus_conclusion.txt","w")
    
    new_config.write("# consensus of the following sequences:\n# ")            # create summary of analysis in own file
    for name in files:
        new_config.write(name.split(".")[0])
        new_config.write(" ")
    new_config.write("\n\n")
    new_config.write("# To modify this configuration, you can add the following line to the blueprint.txt in sRNA_create and follow the instructions there:\n")
    new_config.write("consensus_config;%d" %amount_sequences)
    new_config.write(" |M; %d" %median_richUA)
    new_config.write("; %.2f" %A_percentage)
    new_config.write("; %.2f" %((1-(A_percentage+U_percentage))/2))
    new_config.write("; %.2f" %U_percentage)
    new_config.write(" |G; %d" %median_gap)
    new_config.write(" |S; %d" %int(((min_stemL - free_bp)/2)))
    new_config.write("; %d " %free_bp)
    new_config.write("|M; %d" %median_polyU)
    new_config.write("; %.2f" %((1-(U_percentage_polyU))/3))
    new_config.write("; %.2f" %((1-(U_percentage_polyU))/3))
    new_config.write("; %.2f" %U_percentage_polyU)

    new_config.write("\n\n")
    new_config.write("#######################################\n")
    new_config.write("# amount\n")
    new_config.write("amount = %d\n" %amount_sequences)
    new_config.write("#######################################\n") 
    new_config.write("# UA-rich motif\n")
    new_config.write("#---------------------\n")
    new_config.write("# length\n")
    new_config.write("length_richUA = %d\n\n" %median_richUA)
    new_config.write("# share U on all bases\n")
    new_config.write("U_percentage_richUA = %.2f\n\n" %U_percentage)
    new_config.write("# share A on all bases\n")
    new_config.write("A_percentage_richUA = %.2f\n\n" %A_percentage)
    new_config.write("#---------------------\n")
    new_config.write("# gap motif (distance between stem and UA-motif)\n")
    new_config.write("#---------------------\n")
    new_config.write("# length\n")
    new_config.write("length_gap = %d\n\n" %median_gap)
    new_config.write("#---------------------\n")
    new_config.write("# stemloop motif\n")
    new_config.write("#---------------------\n")
    new_config.write("# length stem (amount basepairs)\n")
    new_config.write("length_stem = %d\n" %int(((min_stemL - free_bp)/2)))
    new_config.write("# length loop (amount unpaired bases)\n")
    new_config.write("length_loop = %d\n" %free_bp)
    new_config.write("#---------------------\n")
    new_config.write("# poly-U- motif\n")
    new_config.write("#---------------------\n")
    new_config.write("# length\n")
    new_config.write("length_polyU = %d\n\n" %median_polyU)
    new_config.write("# share U on all bases\n")
    new_config.write("U_percentage_polyU = %.2f\n\n" %U_percentage_polyU)
    
    new_config.close()
    
    os.chdir('..')
    os.chdir('..')
    path = os.getcwd()                                                         # change working directory
    os.chdir(path +'/sRNAgen')
    
    new_config = open("config_single.yml","w")                                 # overwrite config_single
    
    new_config.write("# set parameters for sequences to be generated\n")
    new_config.write("# consensus of the following sequences:\n# ")
    for name in files:
        new_config.write(name.split(".")[0])
        new_config.write(" ")
    new_config.write("\n\n")
    new_config.write("# To add this configuration to your queue of multiple outputs, add the following line to the config file 'config_multiple.yml':\n")
    new_config.write("consensus_config | %d " %amount_sequences)
    new_config.write("| %d" %median_polyU)
    new_config.write("; %.2f" %U_percentage_polyU)
    new_config.write(" | %d" %median_richUA)
    new_config.write("; %.2f" %U_percentage)
    new_config.write("; %.2f" %A_percentage)
    new_config.write(" | %d" %median_gap)
    new_config.write(" | %d" %int(((min_stemL - free_bp)/2)))
    new_config.write("; %d | " %free_bp)
    for name in files:
        new_config.write(name.split(".")[0])
        new_config.write(" ")
    new_config.write("\n\n")
    new_config.write("#######################################\n")
    new_config.write("# amount\n")
    new_config.write("amount = %d\n" %amount_sequences)
    new_config.write("#######################################\n")
    new_config.write("# UA-rich motif\n")
    new_config.write("#---------------------\n")
    new_config.write("# length\n")
    new_config.write("length_richUA = %d\n\n" %median_richUA)
    new_config.write("# share U on all bases\n")
    new_config.write("U_percentage_richUA = %.2f\n\n" %U_percentage)
    new_config.write("# share A on all bases\n")
    new_config.write("A_percentage_richUA = %.2f\n\n" %A_percentage)
    new_config.write("#---------------------\n")
    new_config.write("# gap motif (distance between stem and UA-motif)\n")
    new_config.write("#---------------------\n")
    new_config.write("# length\n")
    new_config.write("length_gap = %d\n\n" %median_gap)
    new_config.write("#---------------------\n")
    new_config.write("# stemloop motif\n")
    new_config.write("#---------------------\n")
    new_config.write("# length stem (amount basepairs)\n")
    new_config.write("length_stem = %d\n" %int(((min_stemL - free_bp)/2)))
    new_config.write("# length loop (amount unpaired bases)\n")
    new_config.write("length_loop = %d\n" %free_bp)
    new_config.write("#---------------------\n")
    new_config.write("# poly-U- motif\n")
    new_config.write("#---------------------\n")
    new_config.write("# length\n")
    new_config.write("length_polyU = %d\n\n" %median_polyU)
    new_config.write("# share U on all bases\n")
    new_config.write("U_percentage_polyU = %.2f\n\n" %U_percentage_polyU)
    
    new_config.close()
#"""
