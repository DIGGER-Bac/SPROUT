# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 11:40:48 2024

@author: Benedikt Wiebach
Email: bndkt.wbch@gmail.com
"""
import os
import subprocess
from openpyxl import Workbook  
import sys

###############################################################################
# read in configurations from config_multiple
###############################################################################

configuration_file = open("config_multiple.yml","r")                                # read out configuration from file
configuration = configuration_file.readlines()
configuration_file.close()
configurations = []
config_complete = []
try:
    s = -1
    for line in configuration:
        if "#" in line:
            continue
        elif not line.strip():
            continue
        else:
            s = s + 1
            config = line.split("|")
            config_A = config[2].split(";")
            config_B = config[3].split(";")
            config_C = config[5].split(";")
            configurations.append([config[0], int(config[1].strip()), int(config_A[0].strip()), float(config_A[1].strip()), int(config_B[0].strip()), float(config_B[1].strip()), float(config_B[2].strip()), int(config[4].strip()), int(config_C[0].strip()), int(config_C[1].strip()), config[6]])
            config_complete.append(line)
except (ValueError,IndexError):
    print("Configuration #"+str(s+1)+" does not match the expected form. Please change it in 'config_multiple.yml' and try again.")
    sys.exit()
    
###############################################################################
# create sequences from every configuration
###############################################################################

output_sequences = []

for entry in configurations:                                                   # update configuration file
    new_config = open("config_single.yml","w")
    
    
    new_config.write("# set parameters for sequences to be generated\n")
    new_config.write("# consensus of the following sequences:\n# ")

    new_config.write(entry[10])
    new_config.write("\n\n")
    new_config.write("# To add this configuration to your queue of multiple outputs, add the following line to the config file 'config_multiple.yml':\n")
    new_config.write(entry[0])
    new_config.write("| %d " %entry[1])
    new_config.write("| %d" %entry[2])
    new_config.write("; %.2f" %entry[3])
    new_config.write(" | %d" %entry[4])
    new_config.write("; %.2f" %entry[5])
    new_config.write("; %.2f" %entry[6])
    new_config.write(" | %d" %entry[7])
    new_config.write(" | %d" %entry[8])
    new_config.write("; %d | " %entry[9])
    new_config.write(entry[10])
    new_config.write("\n\n")
    new_config.write("#######################################\n")
    new_config.write("# amount\n")
    new_config.write("amount = %d\n" %entry[1])
    new_config.write("#######################################\n")
    new_config.write("# UA-rich motif\n")
    new_config.write("#---------------------\n")
    new_config.write("# length\n")
    new_config.write("length_richUA = %d\n\n" %entry[4])
    new_config.write("# share U on all bases\n")
    new_config.write("U_percentage_richUA = %.2f\n\n" %entry[5])
    new_config.write("# share A on all bases\n")
    new_config.write("A_percentage_richUA = %.2f\n\n" %entry[6])
    new_config.write("#---------------------\n")
    new_config.write("# gap motif (distance between stem and UA-motif)\n")
    new_config.write("#---------------------\n")
    new_config.write("# length\n")
    new_config.write("length_gap = %d\n\n" %entry[7])
    new_config.write("#---------------------\n")
    new_config.write("# stemloop motif\n")
    new_config.write("#---------------------\n")
    new_config.write("# length stem (amount basepairs)\n")
    new_config.write("length_stem = %d\n" %entry[8])
    new_config.write("# length loop  (amount unpaired bases)\n")
    new_config.write("length_loop = %d\n" %entry[9])
    new_config.write("#---------------------\n")
    new_config.write("# poly-U- motif\n")
    new_config.write("#---------------------\n")
    new_config.write("# length\n")
    new_config.write("length_polyU = %d\n\n" %entry[2])
    new_config.write("# share U on all bases\n")
    new_config.write("U_percentage_polyU = %.2f\n\n" %entry[3])
    
    new_config.close()



#     start sRNA_generate
###############################################################################

    path = os.getcwd()                                                         # change working directory
    os.chdir(path +'/utils')
    subprocess.run("python3 sRNA_generate_I.py", shell=True )
    os.chdir('..')

#     read in new generated sequences
###############################################################################

    path = os.getcwd()                                                         # change working directory
    os.chdir(path +'/Output')
    output_file = open("generated_sequences.fasta","r")                                # read out configuration from file
    new_sequences = output_file.readlines()
    output_file.close()
    
    output_sequences.append([entry[0], new_sequences])
    os.chdir('..')

###############################################################################
# organize data for Excel sheet
###############################################################################
path = os.getcwd()                                                             # change working directory
os.chdir(path +'/Output')

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
        tmp_str = set_seq[0].strip() + "-" + str(i+1)
        data_excel.append([tmp_str, tmp_sequences[i]])
    
output_file = open("generated_sequences.fasta","w")                                # read out configuration from file
p = -1
for con in output_sequences:
    p = p + 1
    q = -1
    for line in con:
        q = q + 1
        r = -1
        for seque in line:
            r = r + 1
            if ">" in seque:
                output_file.write(">"+str(p+1)+"_"+output_sequences[p][0].strip()+"_"+output_sequences[p][q][r].split(">")[1].strip()+"\n")
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
for sequence in data_excel:                                                 # append each found occurence ordered by motif
    i = i + 1
    ws.append(sequence)
    
ws.column_dimensions["A"].width = 20                                           # set width of cells so that all information is shown
ws.column_dimensions["B"].width = 60

sheet2 = wb.create_sheet("source_motives")                                     # save information about which motives was searched for
wb.active = wb['source_motives']
sheet2.column_dimensions["A"].width = 250

data_info = ["name of configuration | amount of sequences to be created | length_polyU; U_percentage_polyU | length_richUA; U_percentage_richUA; A_percentage_richUA | length_gap | length_stem (amount basepairs); length_loop | additional data"]
sheet2.append(data_info) 
data_info = ["                                      general data                                                  |                    poly-U motif                        |                                              rich-UA motif                                          |  gap motif  |           stemloop motif         | sRNA names of the consensus"] 
sheet2.append(data_info)   

k = 3
for entry in config_complete:
    k = k + 1
    sheet2.insert_rows(k)
    sheet2["A"+str(k)] = entry
  
wb.active = wb['new_sequences']                                                # set first page as default when opening

while True:
    try:
        wb.save(filename="generated_sequences.xlsx")                           # save file
        break
    except PermissionError:
        print("Can not save the new sequences, because the file 'generated_sequences.xlsx' is open at the moment!")
        print("Please close the file and try again with ENTER.")
        input("WARNING: The file will be overwritten in the process!")
#'''
