# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 12:57:36 2024

@author: Benedikt Wiebach
Email: bndkt.wbch@gmail.com
"""
 
import os
import subprocess
from openpyxl import Workbook  

path = os.getcwd()                                                             # change working directory
os.chdir(path +'/utils')
subprocess.run("python3 sRNA_generate_II.py", shell=True )
os.chdir('..')

output_sequences = []
config_complete =[]
#     read in new generated sequences

config_file = open("config_single.yml","r")   
config = config_file.readlines()
config_file.close()

z = -1
for line in config:
    z = z + 1
    if "|" in line:
        config_complete.append(config[z])

path = os.getcwd()                                                             # change working directory
os.chdir(path +'/Output')

output_file = open("generated_sequences.fasta","r")                                # read out configuration from file
new_sequences = output_file.readlines()
output_file.close()
    
output_sequences.append([config_complete[0].split("|")[0], new_sequences])
os.chdir('..')
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
        data_excel.append([tmp_str, tmp_sequences[i].split(",")[1], int(tmp_sequences[i].split(",")[2].split(":")[1].strip()), tmp_sequences[i].split(",")[0]])

###############################################################################
# save sequences in Excel-sheet
###############################################################################
#'''
path = os.getcwd()                                                             # change working directory
os.chdir(path +'/Output')


wb = Workbook()
ws = wb.active                                                                 # set up workbook
ws.title = "new_sequences"                                                     # create page

data_general = ["ID","ARN motif","gap","sequence"]
ws.append(data_general)                                                        # create headings

i = -1
for sequence in data_excel:                                                 # append each found occurence ordered by motive
    i = i + 1
    ws.append(sequence)
    
ws.column_dimensions["A"].width = 20                                           # set width of cells so that all information is shown
ws.column_dimensions["B"].width = 20
ws.column_dimensions["C"].width =  5
ws.column_dimensions["D"].width = 60

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