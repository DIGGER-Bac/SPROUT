# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 13:15:03 2024

@author: Benedikt Wiebach
Email: bndkt.wbch@gmail.com
"""

##############################################################################

# imports
import subprocess
import os
import sys
from openpyxl import Workbook
from openpyxl.utils import get_column_letter                                 
from openpyxl.formatting.rule import ColorScaleRule, CellIsRule, FormulaRule, Rule

from openpyxl.styles.differential import DifferentialStyle
from openpyxl.styles import Color, PatternFill
import time
import shutil

color_scale_rule = ColorScaleRule(start_type='num', start_value=1, start_color='c65911',
                          mid_type='percent', mid_value=50, mid_color='70ad47',
                          end_type='num', end_value=2, end_color='4472c4')

redFill = PatternFill(start_color='f9cab1',
                   end_color='f9cab1',
                   fill_type='solid')
greenFill = PatternFill(start_color='c7dfb7',
                   end_color='c7dfb7',
                   fill_type='solid')
yellowFill = PatternFill(start_color='ffe694',
                   end_color='ffe694',
                   fill_type='solid')

path = os.getcwd()
try:
    os.mkdir(path+"/Input/processing")
except OSError:
    print("")                                                             # change working directory
#copy perl-script in data folder
shutil.copy("/usr/local/share/ViennaRNA/bin/relplot.pl",path+"/Input/processing/relplot.pl") # copy script of ViennaRNA for positional entropy tothe right position
#shutil.copy("relplot.pl",path+"/Input/processing/relplot.pl")
os.chdir(path +'/Input')                                                        # change working directory

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
            if "," in temp[z+1]:
                new_filename.write(temp[z+1].split(",")[0])
            else:
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
                print("There are no sequences in .fasta format to be analyzed in the folder "+path)
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
for file in files:
    os.remove(file.strip())
os.remove("files.txt")
time.sleep(1)
os.remove("RNAfold.txt")




longest_seq = 0

for sequence in Seq:
    length_sequence = len(sequence)
    if length_sequence > longest_seq:
        longest_seq = length_sequence
    

align_right = longest_seq + 5

###############################################################################
# organize data for Excel sheet
###############################################################################

data_excel_seq = []
data_excel_struct = []
data_excel_prob = []

j = -1
for set_seq in Seq:
    j = j + 1
    data_excel_seq.append([])
    data_excel_struct.append([])
    data_excel_prob.append([])
    i = -1 
    for aligning in range((align_right-len(set_seq.strip()))):
        i = i + 1
        data_excel_seq[j].append(" ")
        data_excel_struct[j].append(" ")
        data_excel_prob[j].append(" ")
    length = -1
    for base in Seq[j].strip():    
        length = length + 1
        data_excel_seq[j].append(base)
    z = -1
    for base in Seq_2DStruct[j]:
        z = z + 1
        if z > length:
            continue
        data_excel_struct[j].append(base)
    for base in Seq_Prob[j]:
        data_excel_prob[j].append(base)
    data_excel_seq[j].append(files[j].split(".")[0])


###############################################################################
# save sequences in Excel-sheet
###############################################################################
#'''
os.chdir('..')
path = os.getcwd() 
shutil.rmtree(path+"/processing")
os.chdir('..')
path = os.getcwd()                                                             # change working directory
os.chdir(path +'/Output')
wb = Workbook()
ws = wb.active                                                                 # set up workbook
ws.title = "new_sequences"                                                     # create page

#data_general = ["ID","sequence"]
#ws.append(data_general)                                                        # create headings

i = -1
for sequence in data_excel_seq:                                                 # append each found occurence ordered by motif
    i = i + 1
    ws.append(list(""))
    ws.append(data_excel_seq[i])
    ws.append(data_excel_struct[i])
    ws.append(list(data_excel_prob[i]))
    ws.append(list(""))
    
    
column = 1    
while column < (align_right+1):
    i = get_column_letter(column)
    ws.column_dimensions[i].width = 2
    column += 1
    



ws.conditional_formatting.add("A1:"+i+"1000", color_scale_rule) # add color rule for each line
#ws.conditional_formatting.add("C"+str(j+2)+":N"+str(j+2), color_scale_rule)


dxf = DifferentialStyle(fill=redFill)
rule = Rule(type="containsText", operator="containsText", text="U", dxf=dxf)
rule.formula = ['NOT(ISERROR(SEARCH("U",A1)))']
ws.conditional_formatting.add("A1:"+i+"1000", rule)

dxf = DifferentialStyle(fill=greenFill)
rule = Rule(type="containsText", operator="containsText", text="A", dxf=dxf)
rule.formula = ['NOT(ISERROR(SEARCH("A",A1)))']
ws.conditional_formatting.add("A1:"+i+"1000", rule)

dxf = DifferentialStyle(fill=yellowFill)
rule = Rule(type="containsText", operator="containsText", text="G", dxf=dxf)
rule.formula = ['NOT(ISERROR(SEARCH("G",A1)))']
ws.conditional_formatting.add("A1:"+i+"1000", rule)





while True:
    try:
        wb.save(filename="visualized_sequences.xlsx")                           # save file
        break
    except PermissionError:
        print("Can not save the new sequences, because the file 'generated_sequences.xlsx' is open at the moment!")
        print("Please close the file and try again with ENTER.")
        input("WARNING: The file will be overwritten in the process!")
#'''

#"""






























