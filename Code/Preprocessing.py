# -*- coding: utf-8 -*-
"""
Created on Thu May 28 16:50:20 2020

@author: minglwu
"""
import pandas as pd
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio import Align

from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment
from collections import defaultdict
from Bio.Alphabet import generic_dna




def match_parent_prodigy(parent_file,prodigy_files):
    parent_dataframe = pd.read_csv('Parents.csv')
    prodigy_dataframe = pd.read_excel('Prodigy.xlsx')
    final_df = pd.merge(prodigy_dataframe, parent_dataframe,left_on='Parent construct',right_on='Name')
    return final_df



def turn_df_into_fasta(merged_df,fasta_output):
    for parent_construct,paired_df in merged_df.groupby('Parent construct'):
        
        starting_seg ='ACCATG'  #GCC
        #print(parent_construct)
        paired_df = paired_df.reset_index(drop= True)
        end_seg  = str(Seq(paired_df.loc[0,'FullSeq'][0:12]).reverse_complement())
        
        parent = paired_df.loc[0,'DNA']
        starting_index = parent.find(starting_seg)
        ending_index = parent.find(end_seg)+12
        parent = parent[starting_index:ending_index]
        
        seq = []
#        sequences = list(paired_df.FullSeq)
        
        seq.append(SeqRecord(Seq(parent),id ='Mother'))
        for i in range(len(paired_df)) :
            seq.append(SeqRecord(Seq(paired_df.loc[i,'FullSeq']).reverse_complement()
                                 , id = paired_df.loc[i,'NGSid']))
            
        SeqIO.write(seq, fasta_output+'/%s.fasta'  %parent_construct, 'fasta')


"""
parent_file ='Parents.csv'
prodigy_files = 'Prodigy.xlsx'

merged_df = match_parent_prodigy(parent_file,prodigy_file)
turn_df_into_fasta(merged_df, fasta_output)
"""