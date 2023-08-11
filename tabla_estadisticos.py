# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 15:50:19 2023

@author: Windows
"""

import pandas as pd
import numpy as np
import obtener_claves

"""Este programa toma como parámetro la matriz libs RPM total y hace la suma d
de la expresión de los genomas, ordenandolos de mayor a menor. También añade
algunos parámetros estadísticos"""

def getStats(inputFile,lim):
    output_dict = {} #las claves son los virus y como valor hay otro diccionario. 
    #Ese diccionario tiene 4 tipos de clave (una para cada tipo de secuencia) y como valor un vector con los distintos valores de media, mediana...
    output_stat = [] #para visualizar la tabla con las estadísticas
    for index,row in inputFile.iterrows():
        #print("going for "+row['gene']+"...")
        if "_" not in str(row['gene']): continue
        ids = row['gene'].split("_",2)
        if ids[0]+"_"+ids[1] not in output_dict.keys():
            output_dict[ids[0]+"_"+ids[1]] = {}
        if len(ids)>2:
            output_dict[ids[0]+"_"+ids[1]][ids[2]] = []
        else:
            ids.append("seq")
            output_dict[ids[0]+"_"+ids[1]]["seq"] = []
        tmp_stat = [] 
        tmp_stat.append(row['gene'])
        values = np.array(row[1:])
        values = np.sort(values[values > lim])
        #calculate mean
        tmp_mean = np.mean(values)
        tmp_stat.append(tmp_mean)
        output_dict[ids[0]+"_"+ids[1]][ids[2]].append(tmp_mean)
        #calculate median
        tmp_median = np.median(values)
        tmp_stat.append(tmp_median)
        output_dict[ids[0]+"_"+ids[1]][ids[2]].append(tmp_median)
        #calculate std
        tmp_std = np.std(values)
        tmp_stat.append(tmp_std)
        output_dict[ids[0]+"_"+ids[1]][ids[2]].append(tmp_std)
        #calculate min
        tmp_min = np.amin(values)
        tmp_stat.append(tmp_min)
        output_dict[ids[0]+"_"+ids[1]][ids[2]].append(tmp_min)
        #calculate max
        tmp_max = np.amax(values)
        tmp_stat.append(tmp_max)
        output_dict[ids[0]+"_"+ids[1]][ids[2]].append(tmp_max)
        #calculate sum
        tmp_sum=np.sum(values)
        tmp_stat.append(tmp_sum)
        output_dict[ids[0]+"_"+ids[1]][ids[2]].append(tmp_sum)
        #calculate p10
        tmp_q10 = np.percentile(values,10)
        tmp_stat.append(tmp_q10)
        output_dict[ids[0]+"_"+ids[1]][ids[2]].append(tmp_q10)
        #calculate p25
        tmp_q25 = np.percentile(values,25)
        tmp_stat.append(tmp_q25)
        output_dict[ids[0]+"_"+ids[1]][ids[2]].append(tmp_q25)
        #calculate p75
        tmp_q75 = np.percentile(values,75)
        tmp_stat.append(tmp_q75)
        output_dict[ids[0]+"_"+ids[1]][ids[2]].append(tmp_q75)
        #calculate p90
        tmp_q90 = np.percentile(values,90)
        tmp_stat.append(tmp_q90)
        output_dict[ids[0]+"_"+ids[1]][ids[2]].append(tmp_q90)
        #calculate p99
        tmp_q99 = np.percentile(values,99)
        tmp_stat.append(tmp_q99)
        output_dict[ids[0]+"_"+ids[1]][ids[2]].append(tmp_q99)
        #number of SRX > limite
        output_dict[ids[0]+"_"+ids[1]][ids[2]].append(len(values))
        tmp_stat.append(len(values))
        #save output
        output_stat.append(tmp_stat) 
    

    output_stat = pd.DataFrame(output_stat,columns=["id","mean","median","std","min","max", "sum", "q10", "q25", "q75", "q90", "q99", "expresion>"+str(lim)])
    unidos = pd.merge(obtener_claves.df, output_stat, on='id', how='inner')
    unidos=unidos.sort_values(by="sum", ascending=False)
    unidos.to_csv("matriz_estadistica_etiqueta_aleat.csv",index=False,sep='\t')
    #output_stat.to_csv("test.out",sep='\t',index=False)
    return output_dict
    return output_stat


# input_path=r"C:\Users\Windows\Desktop\DROSOPHILA\CN_drosophila_resultados\matrix_libs_RPMtotal.tsv"
#input_path=r"C:\Users\Windows\Desktop\matrix_non_human_virus_completo_RPMtotal.tsv"
input_path=r"C:\Users\Windows\Desktop\matrix_human_virus_host_completo_aleat_RPMtotal.tsv"
limite=0

genlib = pd.read_table(input_path,sep='\t')
#genlib_long = pd.melt(genlib, id_vars="gene")
dictStats = getStats(genlib,limite)


