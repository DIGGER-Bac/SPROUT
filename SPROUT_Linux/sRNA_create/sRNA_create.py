# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 12:23:58 2024

@author: Benedikt Wiebach
Email: bndkt.wbch@gmail.com
"""
import os
import subprocess
import random
import time
from openpyxl import Workbook  
import sys

###############################################################################
# read in configurations from blueprint.txt
###############################################################################

configuration_file = open("blueprint.txt","r")                                 # read out configuration from file
configuration = configuration_file.readlines()
configuration_file.close()

configurations = []
config_complete = []
for line in configuration:
    if "#" in line:
        continue
    elif not line.strip():
        continue
    else:
        config = line.split("|")
        configurations.append(config)
        config_complete.append(line)

path = os.getcwd()                                                             # change working directory
os.chdir(path +'/Output')

###############################################################################
# create sequences from every configuration
###############################################################################

output_sequences = []
counter = 0
for entry in configurations:
    counter = counter + 1
    output = open("generated_sequences.fasta","w")
    output.close()
    new_sequences = []
    try:
        name = entry[0].split(";")[0]
    except (ValueError,IndexError):
        print("Configuration #"+str(counter)+" does not match the expected form. There is a problem with the parameter 'name'. Please change it in 'config_multiple.yml' and try again.")
        sys.exit()
    try:
        amount_sequences = int(entry[0].split(";")[1])
    except (ValueError,IndexError):
        print("Configuration #"+str(counter)+" does not match the expected form. There is a problem with the parameter 'amount'. Please change it in 'config_multiple.yml' and try again.")
        sys.exit()

    flag_all_sequences_generated = False
    flag_similarity_found = False
    new = -1
    similarity_overflow = 0
    while flag_all_sequences_generated == False:    
        new = new + 1
        
        input_RNAinverse = open("Seq"+str(new+1)+".in","w")
        input_string1 = ""
        input_string2 = ""
        
        if new == amount_sequences-1:
            flag_all_sequences_generated = True
        i = -1    
        for motif in range((len(entry))):
            flag_no_motif_name = False
            i = i +1
            if i > 0:
                flag_no_motif_name = True
            try:
                motif_sep = entry[i].split(";")
            except (ValueError,IndexError):
                print("Configuration #"+str(counter)+" does not match the expected form. Please change it in 'config_multiple.yml' and try again.")
                sys.exit()
            
            try:
                # motif
                if motif_sep[0] == "M":
                    flag_no_motif_name = False
                    try:
                        prob_A = float(motif_sep[2])
                        prob_C = float(motif_sep[3])
                        prob_U = float(motif_sep[4])
                        if (prob_A + prob_C + prob_U) > 1:
                            print("The probabilities of motif 'M' add up to more than 100%. Please change it in 'config_multiple.yml' and try again.")
                            sys.exit()
                        prob_G = 1 - (prob_A + prob_C + prob_U)
                        len_int = int(motif_sep[1])
                        motif_str = random.choices(['A', 'C', 'G', 'U'], [prob_A, prob_C, prob_G, prob_U], k=len_int)
                        for letter in motif_str:                               # create strings for RNAinverse
                            input_string1 = input_string1 + "."
                            input_string2 = input_string2 + letter.lower()
                    except (ValueError,IndexError):
                        print("Configuration #"+str(counter)+" does not match the expected form. There is a problem with the motif 'M'. Please change it in 'config_multiple.yml' and try again.")
                        sys.exit()
                # gap    
                if motif_sep[0] == "G":
                    flag_no_motif_name = False
                    try:
                        len_gap = int(motif_sep[1])
                        for letter in range(int(len_gap)):                     # random initialisation for gap
                            input_string1 = input_string1 + "."
                            input_string2 = input_string2 + "N"
                    except (ValueError,IndexError):
                        print("Configuration #"+str(counter)+" does not match the expected form. There is a problem with the motif 'G'. Please change it in 'config_multiple.yml' and try again.")
                        sys.exit()
                # Stemloop
                if motif_sep[0] == "S":
                    flag_no_motif_name = False
                    try:
                        length_stem = int(motif_sep[1])
                        free_bp = motif_sep[2]
                        
                        for letter in range(int(length_stem)):
                            input_string1 = input_string1 + "("
                            input_string2 = input_string2 + "N"
                        for letter in range(int(free_bp)):
                            input_string1 = input_string1 + "."
                            input_string2 = input_string2 + "N"
                        for letter in range(int(length_stem)):
                            input_string1 = input_string1 + ")"  
                            input_string2 = input_string2 + "N"
                    except (ValueError,IndexError):
                        print("Configuration #"+str(counter)+" does not match the expected form. There is a problem with the motif 'S'. Please change it in 'config_multiple.yml' and try again.")
                        sys.exit()
                # ARN
                if motif_sep[0] == "A":
                    flag_no_motif_name = False
                    try:
                        prob_A_II = float(motif_sep[2])
                        prob_G_II = 1 - float(motif_sep[2])
                        prob_C_II = 0
                        prob_U_II = 0
                        prob_A_III = float(motif_sep[3])
                        prob_U_III = float(motif_sep[5])
                        prob_C_III = float(motif_sep[4])
                        if (prob_A_III + prob_C_III + prob_U_III) > 1:
                            print("The probabilities of  position III in motif 'A' add up to more than 100%. Please change it in 'config_multiple.yml' and try again.")
                            sys.exit()
                        prob_G_III = (1 - (prob_A_III + prob_U_III + prob_C_III))
            
                    
                        # Generate the chain
                        length = int(int(motif_sep[1])/3)
                        ARN_pos_II = random.choices(['A', 'C', 'G', 'U'], [prob_A_II, prob_C_II, prob_G_II, prob_U_II], k=length)
                        ARN_pos_III = random.choices(['A', 'C', 'G', 'U'], [prob_A_III, prob_C_III, prob_G_III, prob_U_III], k=length)
                
                        output_ARN = ""
                        i = -1
                        for triplet in range(length):
                            i = i + 1
                            output_ARN = output_ARN + "A" + ARN_pos_II[i].upper() + ARN_pos_III[i].upper()
                        for letter in output_ARN:
                            input_string1 = input_string1 + "."
                            input_string2 = input_string2 + letter.lower()
                    except (ValueError,IndexError):
                        print("Configuration #"+str(counter)+" does not match the expected form. There is a problem with the motif 'ARN'. Please change it in 'config_multiple.yml' and try again.")
                        sys.exit()
                # Sequence
                if motif_sep[0] == "Seq":
                    flag_no_motif_name = False
                    try:
                        for letter in motif_sep[1].strip():                            # create strings for RNAinverse
                            if letter == " ":
                                print("Configuration #"+str(counter)+" does not match the expected form. There is a gap in the motif 'Seq'. Please change it in 'config_multiple.yml' and try again.")
                                sys.exit()
                            if letter.lower() != "a" and letter.lower() != "c" and letter.lower() != "g" and letter.lower() != "u":
                                print("Configuration #"+str(counter)+" does not match the expected form. There is a wrong letter in the motif 'Seq'. Please change it in 'config_multiple.yml' and try again.")
                                sys.exit()
                            input_string1 = input_string1 + "."
                            input_string2 = input_string2 + letter.lower()
                    except (ValueError,IndexError):
                        print("Configuration #"+str(counter)+" does not match the expected form. There is a problem with the motif 'Seq'. Please change it in 'config_multiple.yml' and try again.")
                        sys.exit()
                #else
                if flag_no_motif_name == True:
                    print("At least one motif in configuration #"+str(counter)+" is not correctly defined. Please change it in 'config_multiple.yml' and try again.")
                    sys.exit()
            except (ValueError,IndexError):
                print("Configuration #"+str(counter)+" does not match the expected form. There is a problem with the motif naming. Please change it in 'config_multiple.yml' and try again.")
                sys.exit()

        
        input_RNAinverse.write(input_string1)
        input_RNAinverse.write('\n'+input_string2)  
        
        input_RNAinverse.close()
      
        out_RNAinverse = open("RNAinverse.txt","w")
        subprocess.run("RNAinverse < Seq"+str(new+1)+".in -R1 -Fp", stdout=out_RNAinverse, shell=True )
        out_RNAinverse.close()
    
        # add new sequence to data
        #######################################################################
        
        out_RNAinverse = open("RNAinverse.txt","r")                            # read out predicted structure
        output = open("generated_sequences.fasta","a")
        name_sequence = "> " + "sequence_" + str(new+1)
        tmp7 = out_RNAinverse.readlines()
        if tmp7[0] != "":
            new_sequence = tmp7[0].split(" ")
            for seq in new_sequences:                                          # only add if no similarity to other new sequences was found
                if new_sequence[0].upper() == seq:
                    flag_all_sequences_generated = False
                    new = new - 1
                    flag_similarity_found = True
                    similarity_overflow = similarity_overflow + 1
                    print("similarity found and deleted")
                    break
            if flag_similarity_found == False:   
                new_sequences.append(new_sequence[0].upper())
                output.write(name_sequence + "\n")
                output.write(new_sequence[0].upper() + "\n")
                print("generate Seq ",new+1, " / ", amount_sequences)
            elif flag_similarity_found == True:
                if similarity_overflow > 20:                                       # abort if parameters are too strict
                    print("Too strict parameters to generate different sequences! Please extend the parameters or lower the amount of sequences to be generated and try again!")
                    sys.exit()
                flag_similarity_found = False                                  
        out_RNAinverse.close()
        output.close()
        
    time.sleep(2)
    for new in range(amount_sequences):
        os.remove("Seq"+str(new+1)+".in")
    os.remove("RNAinverse.txt")

    output_file = open("generated_sequences.fasta","r")                        # read out configuration from file
    new_sequences = output_file.readlines()
    output_file.close()
    output_sequences.append([entry[0], new_sequences])
os.remove("generated_sequences.fasta") 
   
###############################################################################
# organize data for Excel sheet
###############################################################################

data_excel = []
j = -1
for set_seq in output_sequences:
    j = j + 1
   
    tmp_sequences = []
    k = -1
    for line in set_seq[1]:
        k = k + 1
        if "#" in line:
            continue
        if ">" in line:
            tmp_sequences.append(set_seq[1][k+1].strip())
    i = -1 
    for sequence in range(len(tmp_sequences)):
        i = i + 1
        tmp_str = set_seq[0].split(";")[0].strip() + "-" + str(i+1)
        data_excel.append([tmp_str, tmp_sequences[i]])
       
output_file = open("created_sequences.fasta","w")                                # read out configuration from file
p = -1
for con in output_sequences:
    p = p + 1
    q = -1
    for line in con:
        q = q + 1
        r = -1
        for seque in line:
            r = r + 1
            if ";" in line:
                continue
            if ">" in seque:
                output_file.write(">"+str(p+1)+"_"+output_sequences[p][0].split(";")[0].strip()+"_"+output_sequences[p][q][r].split(">")[1].strip()+"\n")
            else:
                output_file.write(output_sequences[p][q][r])
output_file.close()
    
###############################################################################
# save sequences in Excel-sheet
###############################################################################
#'''
wb = Workbook()
ws = wb.active                                                                 # set up workbook
ws.title = "new_sequences"                                                     # create page

data_general = ["ID","sequence"]
ws.append(data_general)                                                        # create headings

i = -1
for sequence in data_excel:                                                    # append each found occurence ordered by motif
    i = i + 1
    ws.append(sequence)
    
ws.column_dimensions["A"].width = 20                                           # set width of cells so that all information is shown
ws.column_dimensions["B"].width = 60

sheet2 = wb.create_sheet("source_motives")                                     # save information about which motives was searched for
wb.active = wb['source_motives']
sheet2.column_dimensions["A"].width = 250

data_info = ["for help with the short codes please refer to 'blueprint.txt'"] 
sheet2.append(data_info)   

k = 3
for entry in config_complete:
    k = k + 1
    sheet2.insert_rows(k)
    sheet2["A"+str(k)] = entry
  
wb.active = wb['new_sequences']                                                # set first page as default when opening

while True:
    try:
        wb.save(filename="created_sequences.xlsx")                           # save file
        break
    except PermissionError:
        print("Can not save the new sequences, because the file 'generated_sequences.xlsx' is open at the moment!")
        print("Please close the file and try again with ENTER.")
        input("WARNING: The file will be overwritten in the process!")
#'''    