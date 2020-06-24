# -*- coding: utf-8 -*-
"""
Created on Fri May 22 10:51:09 2020

@author: minglwu
"""

import math 
import pandas as pd



def Tm_calculate(c, conc_primer, conc_salt,conc_mg):
    CG_percent = (c.count('C') + c.count('G'))/len(c)
    s = 0
    h = 0 
    # c = Sequence
    # conc_primer in nM (500)
    # conc_salt in mM (50)
    # conc_mg  in mM (2)
     
    enthalpy_table = { 
            'AA':-7.9, 'AC':-8.4, 'AG':-7.8,  "AT":-7.2,
            "CA":-8.5, "CC":-8.0, "CG":-10.6, "CT":-7.8,
            "GA":-8.2, "GC":-9.8, "GG":-8.0, "GT":-8.4,
            "TA":-7.2, "TC":-8.2, "TG":-8.5, "TT": -7.9,
            }
    entropy_table = {
            'AA':-22.2, 'AC':-22.4, 'AG':-21.0,  "AT":-20.4,
            "CA":-22.7, "CC":-19.9, "CG":-27.2, "CT":-21.0,
            "GA":-22.2, "GC":-24.4, "GG":-19.9, "GT":-22.4,
            "TA":-21.3, "TC":-22.2, "TG":-22.7, "TT": -22.2,
            }
    
    salt_effect= (conc_salt/1000)+((conc_mg/1000) * 140)
    s = s + 0.368 * (len(c)-1)* math.log(salt_effect)
    
    if c[0] =="G" or c[0] == "C":
        h = h + 0.1 
        s = s + (-2.8)
    elif  c[0] =="A" or c[0] == "T":
        h = h + 2.3
        s = s + (4.1)
        
    if c[-1] =="G" or c[-1] == "C":   
        h = h + 0.1 
        s = s + (-2.8)
    elif  c[-1] =="A" or c[-1] == "T":
        h = h + 2.3
        s = s + (4.1)
    
    for i in range(len(c)-1):
        #print(c[i:i+2])
        h = h + enthalpy_table[c[i:i+2]]
        s = s + entropy_table[c[i:i+2]]
        
    tm=((1000*h)/(s+(1.987*math.log(conc_primer/2000000000))))-273.15
    return CG_percent, tm



"""
- Annealing region ~15-30 bp 
- Tm ~ 60-70 °C, ∆Tm ≤ 5 °C 
- %GC ~40-60% (can be 0-7% less than 40% and up to ~75% for PCR automation purposes (Not important))
"""
def optimize_best_primer(P1, P2):
    #P1 is forward primer
    #P2 is backward primer 
    
    if P1 == 'No Error':
        # return nothing
        item  = pd.Series({'number of bases': 'N/A','Seq':'N/A',
                                            'CG_percent':'N/A','Tm':'N/A'})
        return item, item 
    
    
    P1_dataframe = pd.DataFrame(columns=['number of bases','Seq','CG_percent','Tm'])
    P2_dataframe = pd.DataFrame(columns=['number of bases','Seq','CG_percent','Tm'])
    
    for i in range(15,len(P1)):
        CG_percent, Tm = Tm_calculate(P1[0:i],500,50,2)
        P1_dataframe = P1_dataframe.append({'number of bases': len(P1[0:i]),'Seq':P1[0:i],
                                            'CG_percent':CG_percent,'Tm':Tm},ignore_index=True)
        CG_percent, Tm = Tm_calculate(P2[0:i],500,50,2)
        P2_dataframe = P2_dataframe.append({'number of bases': len(P2[0:i]),'Seq':P2[0:i],
                                            'CG_percent':CG_percent,'Tm':Tm},ignore_index=True)
        
    
    for i in range(len(P1_dataframe)-1,0-1,-1):
        T1_ok = False 
        T2_ok = False
        deltaT_ok = False
        if (P1_dataframe.loc[i, 'Tm']>60) and (P1_dataframe.loc[i,'Tm']<=70):
            T1_ok = True
            for j in range(len(P2_dataframe)-1,0-1,-1): 
                if abs(P2_dataframe.loc[j, 'Tm'] - P1_dataframe.loc[i, 'Tm'])<5:
                    deltaT_ok = True 
                if P2_dataframe.loc[j,'Tm'] > 60 and P2_dataframe.loc[j,'Tm'] <70:
                    T2_ok = True 
                if  (T1_ok+T2_ok+deltaT_ok)==3  :
                    P1_index = i 
                    P2_index = j
                    break 
            if (T1_ok+T2_ok+deltaT_ok)==3 :
                break
        else:
            continue
        
    if (T1_ok * T2_ok * deltaT_ok) == 0:
        # return nothing
        item  = pd.Series({'number of bases': 'N/A','Seq':'N/A',
                                            'CG_percent':'N/A','Tm':'N/A'})
        return item, item 
    else:
        return P1_dataframe.iloc[P1_index], P2_dataframe.iloc[P2_index]
    
            #check if Tm is 
#optimize_best_primer(P1, P2)
