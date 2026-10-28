# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 13:35:25 2023

@author: Benedikt Wiebach
Email: bndkt.wbch@gmail.com
"""
###############################################################################
import os
import random
import subprocess
import time
import sys

###############################################################################
# read out configuration from file
###############################################################################

os.chdir('..')
configuration_file = open("config_single.yml","r")                                
configuration = configuration_file.readlines()
configuration_file.close

for line in configuration:
    
    if "#" in line:
        continue
    
    # polyU
    if "amount" in line:
        var = line.split("=")
        if var[1].strip().isdigit():
            amount_sequences = int(var[1].strip())
        else:
            print("'amount_sequences' is not defined as an integer in the configuration file. Please change it and start the program again.")
            sys.exit()
        
    if "length_polyU" in line:
        var = line.split("=")
        if var[1].strip().isdigit():
            length_polyU = int(var[1].strip())
        else:
            print("'length_polyU' is not defined as an integer in the configuration file. Please change it and start the program again.")
            sys.exit()

    if "U_percentage_polyU" in line:
        var = line.split("=")
        try:
            U_percentage_polyU = float(var[1].strip())  
        except (ValueError,IndexError):
            print("'U_percentage_polyU' is not defined as a float in the configuration file. Please change it and start the program again.")
            sys.exit()
        
    # ARN
    
    if "length_ARN " in line:
        var = line.split("=")
        if var[1].strip().isdigit():
            length_ARN  = int(var[1].strip())   
        else:
            print("'length_ARN' is not defined as an integer in the configuration file. Please change it and start the program again.")
            sys.exit()
        
    if "A_percentage_II " in line:
        var = line.split("=")
        try:
            A_percentage_II  = float(var[1].strip())
        except (ValueError,IndexError):
            print("'A_percentage_II' is not defined as a float in the configuration file. Please change it and start the program again.")
            sys.exit()        
        
    if "A_percentage_III " in line:
        var = line.split("=")
        try:
            A_percentage_III  = float(var[1].strip())
        except (ValueError,IndexError):
            print("'A_percentage_III' is not defined as a float in the configuration file. Please change it and start the program again.")
            sys.exit()
        
    if "G_percentage_III " in line:
        var = line.split("=")
        try:
            G_percentage_III  = float(var[1].strip())
        except (ValueError,IndexError):
            print("'G_percentage_III' is not defined as a float in the configuration file. Please change it and start the program again.")
            sys.exit()
        
    if "C_percentage_III " in line:
        var = line.split("=")
        try:
            C_percentage_III  = float(var[1].strip())
        except (ValueError,IndexError):
            print("'C_percentage_III' is not defined as a float in the configuration file. Please change it and start the program again.")
            sys.exit()
    
    # gap    
    if "length_gap" in line:
        var = line.split("=")
        if var[1].strip().isdigit():
            length_gap = int(var[1].strip())   
        else:
            print("'length_gap' is not defined as an integer in the configuration file. Please change it and start the program again.")
            sys.exit()
        
    # stemL 
    if "length_stem" in line:
        var = line.split("=")
        if var[1].strip().isdigit():
            length_stem = int(var[1].strip())   
        else:
            print("'length_stem' is not defined as an integer in the configuration file. Please change it and start the program again.")
            sys.exit()
        
    if "length_loop" in line:
        var = line.split("=")
        if var[1].strip().isdigit():
            free_bp = int(var[1].strip())   
        else:
            print("'length_loop' is not defined as an integer in the configuration file. Please change it and start the program again.")
            sys.exit()

if (A_percentage_III + G_percentage_III + C_percentage_III) > 1:
    print("The percentages for the ARN motif add up to more than 100%! Please change it and start the program again.")
    sys.exit()
length_stemL = length_stem + free_bp
path = os.getcwd()                                                             # change working directory
os.chdir(path +'/Output')

###############################################################################
# write Config to Output
###############################################################################

output = open("generated_sequences.fasta","w")
output.write("# motive length (f.l.t.r.): %d, " %length_ARN)
output.write("%d, " %length_gap)
output.write("%d, " %length_stemL)
output.write("%d" %length_polyU)
output.write(", possibilities: A_II: %.2f" %A_percentage_II)
output.write(", A_III:%.2f" %A_percentage_III)
output.write(", G_III:%.2f" %G_percentage_III)
output.write(", C_III:%.2f" %C_percentage_III)
output.write(", U(polyU): %.2f" %U_percentage_polyU)
output.write(", length loop: %d\n" %free_bp)

print("generate new sequences")
output.close()

# create datastructure to check for similarity
new_sequences = []
new_sequences_ARN = []
flag_all_sequences_generated = False
flag_similarity_found = False
new = -1
similarity_overflow = 0

while flag_all_sequences_generated == False:    
    new = new + 1
    if new == amount_sequences-1:
        flag_all_sequences_generated = True
    
# poly-U
##############################################################################

    prob_A = prob_C = prob_G = (1 - U_percentage_polyU)/3
    prob_U = U_percentage_polyU
        
    len_int = int(length_polyU)
    polyU = random.choices(['A', 'C', 'G', 'U'], [prob_A, prob_C, prob_G, prob_U], k=len_int)
  
# ARN
##############################################################################

    prob_A_II = A_percentage_II
    prob_G_II = 1 - A_percentage_II
    prob_C_II = 0
    prob_U_II = 0
    prob_A_III = A_percentage_III
    prob_G_III = G_percentage_III
    prob_C_III = C_percentage_III
    prob_U_III = (1 - (A_percentage_III + G_percentage_III + C_percentage_III))

    length = int(length_ARN/3)
    ARN_pos_II = random.choices(['A', 'C', 'G', 'U'], [prob_A_II, prob_C_II, prob_G_II, prob_U_II], k=length)
    ARN_pos_III = random.choices(['A', 'C', 'G', 'U'], [prob_A_III, prob_C_III, prob_G_III, prob_U_III], k=length)
    
    output_ARN = ""
    i = -1
    for triplet in range(length):
        i = i + 1
        output_ARN = output_ARN + "A" + ARN_pos_II[i].upper() + ARN_pos_III[i].upper()

# stemloop and concatenation
##############################################################################

    input_RNAinverse = open("Seq"+str(new+1)+".in","w")
    input_string1 = ""
    input_string2 = ""
#    for letter in richUA:                                                      # not needed for Class II
#        input_string1 = input_string1 + "."
#        input_string2 = input_string2 + letter.lower()
#    for letter in range(int(length_gap)):                                      # gap is only added as information
#        input_string1 = input_string1 + "."                                    # -> no need for initialisation
#        input_string2 = input_string2 + "N"
    for letter in range(int((length_stem))):
        input_string1 = input_string1 + "("
        input_string2 = input_string2 + "N"
    for letter in range(int(free_bp)):
        input_string1 = input_string1 + "."
        input_string2 = input_string2 + "N"
    for letter in range(int((length_stem))):
        input_string1 = input_string1 + ")"  
        input_string2 = input_string2 + "N"
    for letter in polyU:
        input_string1 = input_string1 + "."
        input_string2 = input_string2 + letter.lower()
    input_RNAinverse.write(input_string1)
    input_RNAinverse.write('\n'+input_string2)  
    
    input_RNAinverse.close()
    
    out_RNAinverse = open("RNAinverse.txt","w")
    subprocess.run("RNAinverse < Seq"+str(new+1)+".in -R1 -Fp", stdout=out_RNAinverse, shell=True )
    out_RNAinverse.close()

# add new sequence to data
##############################################################################
    
    out_RNAinverse = open("RNAinverse.txt","r")                                # read out predicted structure
    output = open("generated_sequences.fasta","a")
    name_sequence = ">" + "sequence_" + str(new+1)
    tmp7 = out_RNAinverse.readlines()
    if tmp7[0] != "":
        new_sequence = tmp7[0].split(" ")
        for seq in new_sequences:                                              # only add if no similarity to other new sequences was found
            if new_sequence[0].upper() == seq:
                
                flag_all_sequences_generated = False
                new = new - 1
                flag_similarity_found = True
                similarity_overflow = similarity_overflow + 1
                print("similarity found and deleted")
                break
        for seq in new_sequences_ARN:                                              
            if flag_similarity_found == True:
                break
            if output_ARN == seq:
                
                flag_all_sequences_generated = False
                new = new - 1
                flag_similarity_found = True
                similarity_overflow = similarity_overflow + 1
                print("similarity (ARN) found and deleted")
                break
            
        if flag_similarity_found == False:   
            new_sequences.append(new_sequence[0].upper())
            new_sequences_ARN.append(output_ARN)
            output.write(name_sequence + "\n")
            output.write(new_sequence[0].upper() + ", " + output_ARN + ", gap: " + str(length_gap) + "\n")
            print("generate Seq ",new+1, " / ", amount_sequences)
        elif flag_similarity_found == True:
            if similarity_overflow > 20:                                       # abort if parameters are too strict
                print("Too strict parameters to generate different sequences! Please extend the parameters or lower the amount of sequences to be generated and try again!")
                sys.exit()
            flag_similarity_found = False                                      
    out_RNAinverse.close()
    output.close()
    
##############################################################################
# delete unnecessary data
##############################################################################

time.sleep(2)
for new in range(amount_sequences):
    os.remove("Seq"+str(new+1)+".in")
os.remove("RNAinverse.txt")