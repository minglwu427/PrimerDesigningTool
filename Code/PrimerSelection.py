# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio import Align
from Bio.pairwise2 import format_alignment
from collections import defaultdict
import pandas as pd
from Bio.Alphabet import generic_dna
import numpy as np
from Tm_Calculation import *
import os





"""
Open a fastq file and covert it Translation and turn it into output
"""

def Extract(lst): 
    return [item[0] for item in lst] 

table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', } 


inv_map = defaultdict(list)
{inv_map[v].append(k) for k, v in table.items()}






def translation_to_AA(fasta_file, final_dataframe):
    translate = []
    NGS_ID = []
    OG_SEQ = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences = record.translate()
            sequences.id = record.id
            sequences.name = record.name
            translate.append(sequences)
            NGS_ID.append(record.id)
            OG_SEQ.append(str(record.seq))
    final_dataframe['NGS_ID'] = NGS_ID
    final_dataframe['OG SEQ'] = OG_SEQ
    SeqIO.write(translate, "AA_translate.fasta", "fasta")
    return translate, final_dataframe





def align_AA(translate, final_dataframe):
    OG_AA = []
    PP_AA = []
    mother = translate[0]
    for AA_seq in translate:
        pair_result = pairwise2.align.globalxx(str(mother.seq), str(AA_seq.seq))
        #alignments.append([pair_result[0][0], pair_result[0][1]])
        OG_AA.append(pair_result[0][0])
        PP_AA.append(pair_result[0][1])
        
    final_dataframe['OG AA'] = OG_AA
    final_dataframe['PP AA'] = PP_AA
    
    return final_dataframe






def define_mutation(OG_Seq, PP_Seq):
    if OG_Seq.count('-')  == 0 and PP_Seq.count('-')  != 0 :
        return 'Deletion'
    elif OG_Seq.count('-') !=0 and PP_Seq.count('-')==0:
        return 'Insertion'
    elif OG_Seq.count('-') !=0 and PP_Seq.count('-') !=0 and OG_Seq.count('-')<PP_Seq.count('-'):
        return "Deletion and Substitution"
    elif OG_Seq.count('-') !=0 and PP_Seq.count('-') !=0 and OG_Seq.count('-')>PP_Seq.count('-'):
        return 'Insertion and Substitution'
    elif OG_Seq.count('-') !=0 and PP_Seq.count('-') !=0 and OG_Seq.count('-')==PP_Seq.count('-'):
        return  "Substitution"
    elif OG_Seq.count('-')==0 and PP_Seq.count('-')==0 :
        return "No Mutations"
   
def assign_mutation(final_dataframe):
    definition=[]
    for i in range(len(final_dataframe)):
        definition.append(define_mutation(final_dataframe.loc[i, 'OG AA'], 
                                          final_dataframe.loc[i,'PP AA']))
    final_dataframe['Mutation Type'] = definition
    return final_dataframe





def reverse_translation(final_dataframe):
    """
    This script doesn't handle deletion with substitution just yet.  

    """
    OG_translated_DNA = []
    PP_translated_DNA= []
    for i in range(len(final_dataframe)): 
        OG_DNA_RT = ""
        PP_DNA_RT = ""
        OG_nuc = final_dataframe.loc[0, 'OG SEQ']
        PP_nuc = final_dataframe.loc[i,'OG SEQ']
        OG_AA_aligned = final_dataframe.loc[i,'OG AA']
        PP_AA_aligned = final_dataframe.loc[i,'PP AA']
        Mutation = final_dataframe.loc[i,'Mutation Type']
        if Mutation == 'Deletion':
            for j in range(len(OG_AA_aligned)):
                if OG_nuc[3*j: 3*j+3 ] == PP_nuc[3*j: 3*j+3 ]: # handle no mistake
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                elif (OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ]) and OG_AA_aligned[j] == PP_AA_aligned : #handle silent Mutation
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3]
                elif OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ] and PP_AA_aligned[j] == '-' :
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + "---"
                    PP_nuc = PP_nuc[:j] + '---' + PP_nuc[j:]
            OG_translated_DNA.append(OG_DNA_RT)
            PP_translated_DNA.append(PP_DNA_RT)
            
        if Mutation == 'Insertion':
            for j in range(len(PP_AA_aligned)):
                if OG_nuc[3*j: 3*j+3 ] == PP_nuc[3*j: 3*j+3 ]:
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                elif (OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ]) and OG_AA_aligned[j] == PP_AA_aligned[j] : 
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3]
                elif OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ] and OG_AA_aligned[j] == '-' : #handle insertion
                    OG_DNA_RT = OG_DNA_RT + "---"#OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3 ]
                    OG_nuc = OG_nuc[:j] + '---' + OG_nuc[j:]
            OG_translated_DNA.append(OG_DNA_RT)
            PP_translated_DNA.append(PP_DNA_RT)
        if Mutation == 'Substitution':
            OG_AA_aligned = OG_AA_aligned.replace('-','')
            PP_AA_aligned = PP_AA_aligned.replace('-','')
            for j in range(len(PP_AA_aligned)):
                if OG_nuc[3*j: 3*j+3 ] == PP_nuc[3*j: 3*j+3 ]:
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                elif (OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ]) and OG_AA_aligned[j] == PP_AA_aligned[j] : 
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3]
                else: #as long as they are reverse translating its own line everything should be fine
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3]
            OG_translated_DNA.append(OG_DNA_RT)
            PP_translated_DNA.append(PP_DNA_RT)
        if Mutation == 'Insertion and Substitution':
            OG_AA_aligned.count('-')
            OG_AA_aligned = OG_AA_aligned.replace('-','',PP_AA_aligned.count('-'))
            PP_AA_aligned = PP_AA_aligned.replace('-','',PP_AA_aligned.count('-'))
            for j in range(len(PP_AA_aligned)):
                if OG_nuc[3*j: 3*j+3 ] == PP_nuc[3*j: 3*j+3 ]:
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + OG_nuc[3*j: 3*j+3 ]
                elif (OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ]) and OG_AA_aligned[j] == PP_AA_aligned[j] : 
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3]
                elif OG_nuc[3*j: 3*j+3 ] != PP_nuc[3*j: 3*j+3 ] and OG_AA_aligned[j] == '-' : #handle insertion
                    OG_DNA_RT = OG_DNA_RT + "---"#OG_nuc[3*j: 3*j+3 ]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3 ]
                    OG_nuc = OG_nuc[:j] + '---' + OG_nuc[j:]
                else: #as long as they are reverse translating its own line everything should be fine
                    OG_DNA_RT = OG_DNA_RT + OG_nuc[3*j: 3*j+3]
                    PP_DNA_RT = PP_DNA_RT + PP_nuc[3*j: 3*j+3]
            OG_translated_DNA.append(OG_DNA_RT)
            PP_translated_DNA.append(PP_DNA_RT)
        if Mutation == 'No Mutations' or Mutation == "Deletion and Substitution":
            OG_DNA_RT = OG_nuc
            PP_DNA_RT = OG_nuc
            OG_translated_DNA.append(OG_DNA_RT)
            PP_translated_DNA.append(PP_DNA_RT)

    final_dataframe['OG_translated_DNA'] = OG_translated_DNA
    final_dataframe['PP_translated_DNA'] = PP_translated_DNA
    return final_dataframe




def generate_raw_primers(final_dataframe):
    forward_primer_list = []
    backward_primer_list = []
    mutation_region_list = []
    """
    generate_raw_primers
    """
    for i in range(len(final_dataframe)): 
        index = list(final_dataframe.iloc[i])
        first_error_location = 0
        find_first_error = False
        find_last_error = False
        last_error_location = 0
        OG_DNA_RT = final_dataframe.loc[i,'OG_translated_DNA']
        PP_DNA_RT = final_dataframe.loc[i,'PP_translated_DNA']
        if final_dataframe.loc[i,'OG AA'].count('-') > final_dataframe.loc[i,'PP AA'].count('-'):
            dash_num = final_dataframe.loc[i,'PP AA'].count('-')
        else:
            dash_num = final_dataframe.loc[i,'OG_translated_DNA'].count('-')
        OG_AA_aligned = final_dataframe.loc[i,'OG AA'].replace('-','',dash_num)
        PP_AA_aligned = final_dataframe.loc[i,'PP AA'].replace('-','',dash_num)
        for j in range(len(OG_AA_aligned)):
            if find_first_error == True and find_last_error == True:
                break
            elif (OG_AA_aligned[j:j+50] == PP_AA_aligned [j:j+50]) and find_first_error == True and find_last_error == False:
                find_last_error = True
                last_error_location = j
            elif (OG_AA_aligned[j] != PP_AA_aligned [j]) and find_first_error == False: 
                find_first_error = True
                first_error_location = j
            elif (OG_AA_aligned[j] == PP_AA_aligned[j]):
                continue

        first_error_location = first_error_location * 3
        last_error_location = last_error_location * 3
        
        """
        fix edge case
        """
        for i in range(3):
            if PP_DNA_RT[first_error_location+i] == OG_DNA_RT[first_error_location+i]:
                first_error_location = first_error_location +1
            if PP_DNA_RT[last_error_location+i] == OG_DNA_RT[last_error_location+i]:
                first_error_location = first_error_location +1
        """
        assign to the dataframe
        """
        if find_first_error == True and find_last_error == True :
            MM_region = PP_DNA_RT[first_error_location:last_error_location].replace('-',"")
            backward_primer = PP_DNA_RT[first_error_location-30: first_error_location]
            backward_primer = str(Seq(backward_primer).reverse_complement()) #reverse compliment
            forward_primer  = PP_DNA_RT[last_error_location: last_error_location+30]
            mutation_region_list.append(MM_region)
            forward_primer_list.append(forward_primer)
            backward_primer_list.append(backward_primer)
        else:
            mutation_region_list.append('No Error')
            forward_primer_list.append('No Error')
            backward_primer_list.append('No Error')
            
    final_dataframe['forward_raw_primer'] = forward_primer_list
    final_dataframe['backward_raw_primer'] = backward_primer_list
    final_dataframe['mutation_region'] = mutation_region_list
    return final_dataframe
    




def make_best_primer(final_dataframe):
    for i in range(len(final_dataframe)): 
        forward_raw_primer = final_dataframe.loc[i, 'forward_raw_primer']
        backward_raw_primer = final_dataframe.loc[i, 'forward_raw_primer']
        best_forward_primer, best_backward_primer = optimize_best_primer(forward_raw_primer, 
                                                                         backward_raw_primer)
    
        final_dataframe.loc[i,'BFP_number_of_bases'] = best_forward_primer[0]
        final_dataframe.loc[i,'BFP_sequence'] = best_forward_primer[1]
        final_dataframe.loc[i,'BFP_CG_percent'] = best_forward_primer[2]
        final_dataframe.loc[i,'BFP_Tm'] = best_forward_primer[3]
        
        final_dataframe.loc[i,'BBP_number_of_bases'] = best_backward_primer[0]
        final_dataframe.loc[i,'BBP_sequence'] = best_backward_primer[1]
        final_dataframe.loc[i,'BBP_CG_percent'] = best_backward_primer[2]
        final_dataframe.loc[i,'BBP_Tm'] = best_backward_primer[3]
    return final_dataframe    



def adding_overhang_to_primer(final_dataframe):
    for i in range(len(final_dataframe)): 
        if final_dataframe.loc[i ,'mutation_region' ] == 'No Error':
            final_dataframe.loc[i, 'final_forward_primer']  = 'N/A'
            final_dataframe.loc[i, 'final_backward_primer']  = 'N/A'
        else:
            mutation_region = final_dataframe.loc[i ,'mutation_region' ]
            mutation_midpoint = int(len(mutation_region)/2)
            backward_overhang = mutation_region[0:mutation_midpoint]
            backward_overhang = str(Seq(backward_overhang).reverse_complement())
            final_dataframe.loc[i, 'final_backward_primer'] = backward_overhang +final_dataframe.loc[i, 'BFP_sequence']
            
            forward_overhang = mutation_region[mutation_midpoint:]
            final_dataframe.loc[i, 'final_forward_primer']  = forward_overhang + final_dataframe.loc[i, 'BFP_sequence'] 
    return final_dataframe


"""
file_list=os.listdir("C:/Users/minglwu/Desktop/Macro Scenarios/Output/fasta")
for file_name in file_list:
    fasta_file = 'Output/fasta/'+ file_name
    final_dataframe = pd.DataFrame()
    translate, final_dataframe = translation_to_AA(fasta_file,final_dataframe)
    final_dataframe = align_AA(translate,final_dataframe)
    
    final_dataframe = assign_mutation(final_dataframe)
    final_dataframe = reverse_translation(final_dataframe)
    final_dataframe = generate_raw_primers(final_dataframe)
    final_dataframe = make_best_primer(final_dataframe)
    
    final_dataframe = adding_overhang_to_primer(final_dataframe)
    
    
    final_dataframe.to_excel('Output/primer/%s_primer.xlsx' %file_name[0:5])
"""


