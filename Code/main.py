
from Preprocessing import *
#import Tm_Calculation
from PrimerSelection import *
import tkinter as tk
from tkinter import filedialog 
import os 

cwd = os.getcwd()
def getting_input():
    tk.messagebox.showinfo('Message 1', 'Please enter your prodigy file (in Excel)')
    prodigy_file =  filedialog.askopenfilename(initialdir = cwd, title = "Select prodigy file",
                                             filetypes = (('excel files','*.xlsx'),("csv files","*.csv"),("all files","*.*")))
    tk.messagebox.showinfo('Message 2', 'Please enter your parent file (in CSV)')
    parent_file =  filedialog.askopenfilename(initialdir = cwd, title = "Select parent file",
                                             filetypes = (('excel files','*.xlsx'),("csv files","*.csv"),("all files","*.*")))
    
    tk.messagebox.showinfo('Message 3', 'Please enter your fasta output directory')
    fasta_output = filedialog.askdirectory(initialdir = cwd, title  = 'Select fasta output directory')
    
    tk.messagebox.showinfo('Message 4', 'Please enter your final primer output directory')
    primer_output = filedialog.askdirectory(initialdir = cwd, title  = 'Select  final primer output directory')
    
    merged_df = match_parent_prodigy(parent_file,prodigy_file)
    turn_df_into_fasta(merged_df, fasta_output)

    return fasta_output, primer_output



fasta_output, primer_output = getting_input()

file_list=os.listdir(fasta_output)

for file_name in file_list:
    fasta_file = fasta_output + '/' +file_name
    final_dataframe = pd.DataFrame()
    translate, final_dataframe = translation_to_AA(fasta_file,final_dataframe)
    final_dataframe = align_AA(translate,final_dataframe)
    
    final_dataframe = assign_mutation(final_dataframe)
    final_dataframe = reverse_translation(final_dataframe)
    final_dataframe = generate_raw_primers(final_dataframe)
    final_dataframe = make_best_primer(final_dataframe)
    
    final_dataframe = adding_overhang_to_primer(final_dataframe)
    
    
    final_dataframe.to_excel(primer_output + '/%s_primer.xlsx' %file_name[0:5])