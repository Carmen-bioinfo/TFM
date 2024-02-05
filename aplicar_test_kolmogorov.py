# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 12:23:32 2023

@author: Windows
"""
import pandas as pd
import numpy as np

from scipy import stats

#---------------- ESTA PARTE DEL CÓDIGO SOLO LEE Y PREPARA LOS DATOS PARA PODER HACER LA REPRESENTACIÓN
input_path=r"C:\Users\Windows\Desktop\TFM\MIS PROGRAMAS\SRX\tabla_sin_filtrar_hsa.csv"

tabla = pd.read_csv(input_path, header=0, delimiter="\t")

# seleccion=tabla.head(20)

seleccion=tabla.copy()
seleccion = seleccion.rename(columns={'id': 'gene'})
split_genlib=seleccion['gene'].str.split('_', n=2, expand=True)
nombre_NC=split_genlib[0]+"_"+split_genlib[1]
 
tipo_seq=split_genlib[2]

tipo_seq = tipo_seq.replace(to_replace=np.nan, value='seq')

del(seleccion["gene"])

seleccion.loc[:, "tipo_seq"] = tipo_seq
seleccion.loc[:, "gene"] = nombre_NC

# Obtener el índice actual de las columnas
cols = seleccion.columns.tolist()

# Mover la columna "col3" al primer lugar
cols.insert(0, cols.pop(cols.index('gene')))
cols.insert(1, cols.pop(cols.index('tipo_seq')))

# # Reindexar el dataframe con el nuevo orden de columnas
df_trabajo = seleccion.reindex(columns=cols)

#-------APLICAMOS EL TEST PARA COMPARAR ENTRE LAS SECUENCIAS + Y - Y SUS RESPECTIVAS ALEATORIZADAS--------- 
print("-------TEST PARA COMPARAR ENTRE LAS SECUENCIAS + Y - Y SUS RESPECTIVAS ALEATORIZADAS---------") 

plus_seq=df_trabajo[df_trabajo["tipo_seq"]=="seq"] 
plus_seq_sum= np.array(plus_seq['sum'])

random_plus_seq=df_trabajo[df_trabajo["tipo_seq"]=="seq_aleatoria"] 
random_plus_seq_sum=np.array(random_plus_seq['sum'])

statistic, p_value = stats.mannwhitneyu(plus_seq_sum, random_plus_seq_sum)
print("Test para secuencia directa y su aleatorizada")
print("Estadística de prueba:", statistic)
print("Valor p:", p_value)
print("\n")


# --------Selección de la revseq
minus_seq=df_trabajo[df_trabajo["tipo_seq"]=="seq_reverso_complementaria"] 
minus_seq_sum= np.array(minus_seq['sum'])

random_minus_seq=df_trabajo[df_trabajo["tipo_seq"]=="seq_aleatoria_complementaria"] 
random_minus_seq_sum=np.array(random_minus_seq['sum'])

statistic, p_value = stats.mannwhitneyu(minus_seq_sum, random_minus_seq_sum)
print("Test para secuencia complementaria y su aleatorizada")
print("Estadística de prueba:", statistic)
print("Valor p:", p_value)
print("\n")



#------AÑADIENDO VIRUS NO HUMANOS MAPEADOS FRENTE A HUMANOS-------
print("------AÑADIENDO VIRUS NO HUMANOS MAPEADOS FRENTE A HUMANOS-------")

input_path2=r"C:\Users\Windows\Desktop\TFM\MIS PROGRAMAS\SRX\matriz_estadistica_etiqueta_non_human.csv"

tabla_non_human = pd.read_csv(input_path2, header=0, delimiter="\t")
#---primero 
# seleccion=tabla.head(20)

seleccion_non_human=tabla_non_human.copy()

split_genlib_non_human=seleccion_non_human['id'].str.split('_', n=2, expand=True)
nombre_NC_non_human=split_genlib_non_human[0]+"_"+split_genlib_non_human[1]
 
tipo_seq_non_human=split_genlib_non_human[2]

tipo_seq_non_human = tipo_seq_non_human.replace(to_replace=np.nan, value='seq')

del(seleccion_non_human["id"])

seleccion_non_human.loc[:, "tipo_seq"] = tipo_seq_non_human
seleccion_non_human.loc[:, "gene"] = nombre_NC_non_human

# Obtener el índice actual de las columnas
cols = seleccion_non_human.columns.tolist()

# Mover la columna "col3" al primer lugar
cols.insert(0, cols.pop(cols.index('gene')))
cols.insert(1, cols.pop(cols.index('tipo_seq')))

# # Reindexar el dataframe con el nuevo orden de columnas
df_trabajo_non_human = seleccion_non_human.reindex(columns=cols)

#--------Aplicamos test K con directa---------
plus_seq_NH=df_trabajo_non_human[df_trabajo_non_human["tipo_seq"]=="seq"]
plus_seq_NH_sum=np.array(plus_seq_NH['sum'])
statistic, p_value = stats.mannwhitneyu(plus_seq_sum, plus_seq_NH_sum)
print("---Secuencia + y aquella en non human virus---")
print("Estadística de prueba:", statistic)
print("Valor p:", p_value)
print("\n")

#---------Aplicamos test K con revseq----------------
minus_seq_NH=df_trabajo_non_human[df_trabajo_non_human["tipo_seq"]=="seq_reverso_complementaria"]

minus_seq_NH_sum=np.array(minus_seq_NH['sum'])

statistic, p_value = stats.mannwhitneyu(minus_seq_NH_sum, minus_seq_sum)
print("---Secuencia - y aquella en non human virus---")
print("Estadística de prueba:", statistic)
print("Valor p:", p_value)
print("\n")

#-----Test de Kolmogorov entre seq aleatoria directa y seq plus non human virus
statistic, p_value = stats.mannwhitneyu(plus_seq_NH_sum, random_plus_seq_sum)
print("---Entre seq aleatoria directa y seq plus non human virus---")
print("Estadística de prueba:", statistic)
print("Valor p:", p_value)
print("\n")

#-----Test de Kolmogorov entre seq aleatoria minus y seq minus non human virus

statistic, p_value = stats.mannwhitneyu(minus_seq_NH_sum, random_minus_seq_sum)
print("---Entre seq aleatoria minus y seq minus non human virus---")
print("Estadística de prueba:", statistic)
print("Valor p:", p_value)
print("\n")

#--------APLICAMOS EL TEST A LAS DOS ALEATORIZADAS--------
print("--------APLICAMOS EL TEST A LAS DOS ALEATORIZADAS--------")



input_path3=r"C:\Users\Windows\Desktop\TFM\MIS PROGRAMAS\SRX\matriz_estadistica_etiqueta_aleat.csv"

tabla_aleat = pd.read_csv(input_path3, header=0, delimiter="\t")

# seleccion=tabla.head(20)

seleccion_aleat=tabla_aleat.copy()

split_genlib_aleat=seleccion_aleat['id'].str.split('_', n=2, expand=True)
nombre_NC_aleat=split_genlib_aleat[0]+"_"+split_genlib_aleat[1]
 
tipo_seq_aleat=split_genlib_aleat[2]

tipo_seq_aleat = tipo_seq_aleat.replace(to_replace=np.nan, value='seq')

del(seleccion_aleat["id"])

seleccion_aleat.loc[:, "tipo_seq"] = tipo_seq_aleat
seleccion_aleat.loc[:, "gene"] = nombre_NC_aleat

# Obtener el índice actual de las columnas
cols_aleat = seleccion_aleat.columns.tolist()

# Mover la columna "col3" al primer lugar
cols_aleat.insert(0, cols_aleat.pop(cols_aleat.index('gene')))
cols_aleat.insert(1, cols_aleat.pop(cols_aleat.index('tipo_seq')))

# # Reindexar el dataframe con el nuevo orden de columnas
df_trabajo_aleat = seleccion_aleat.reindex(columns=cols_aleat)

random_random_minus=df_trabajo_aleat[df_trabajo_aleat["tipo_seq"]=="seq_aleatoria_complementaria"] 
random_random_minus_sum=np.array(random_random_minus['sum'])
statistic, p_value = stats.mannwhitneyu(random_random_minus_sum, random_minus_seq_sum)
print("---Entre seq aleatoria complementaria y otra seq aleatoria---")
print("Estadística de prueba:", statistic)
print("Valor p:", p_value)
print("\n")


random_random_plus=df_trabajo_aleat[df_trabajo_aleat["tipo_seq"]=="seq_aleatoria"] 
random_random_plus_sum=np.array(random_random_plus['sum'])
statistic, p_value = stats.mannwhitneyu(random_random_plus_sum, random_plus_seq_sum)
print("---Entre seq aleatoria directa y otra seq aleatoria---")
print("Estadística de prueba:", statistic)
print("Valor p:", p_value)
print("\n")

#----APLICAMOS LA PRUEBA CON LAS TABLAS YA FILTRADAS-------
# print("----PRUEBA PARA LAS TABLAS FILTRADAS----")
# input_path4=r"C:\Users\Windows\Desktop\TFM\MIS PROGRAMAS\SRX\tabla_filtrada_hsa.csv"
# tabla_filtrada = pd.read_csv(input_path4, header=0, delimiter="\t")
# f_plus_seq=tabla_filtrada[tabla_filtrada["Tipo_seq"]=="seq"] 
# f_plus_seq_sum= np.array(f_plus_seq['sum'])


print("\n")
print("------T DE STUDENT------")
#----------T DE STUDENT-------
t_statistic, p_value = stats.ttest_ind(plus_seq_sum, random_plus_seq_sum)
print("Test para seq + y su aleatoria")
print("Estadística de prueba (t):", t_statistic)
print("Valor p:", p_value)

t_statistic, p_value = stats.ttest_ind(minus_seq_sum, random_minus_seq_sum)
print("Test para seq - y su aleatoria")
print("Estadística de prueba (t):", t_statistic)
print("Valor p:", p_value)

t_statistic, p_value = stats.ttest_ind(plus_seq_sum, plus_seq_NH_sum)
print("Test para seq + y seq + non human virus")
print("Estadística de prueba (t):", t_statistic)
print("Valor p:", p_value)

t_statistic, p_value = stats.ttest_ind(minus_seq_sum, minus_seq_NH_sum)
print("Test para seq - y seq - non human virus")
print("Estadística de prueba (t):", t_statistic)
print("Valor p:", p_value)

t_statistic, p_value = stats.ttest_ind(plus_seq_NH_sum, random_plus_seq_sum)
print("seq + aleatoria y seq + non human virus")
print("Estadística de prueba (t):", t_statistic)
print("Valor p:", p_value)

t_statistic, p_value = stats.ttest_ind(minus_seq_NH_sum, random_minus_seq_sum)
print("seq - y seq - non human virus")
print("Estadística de prueba (t):", t_statistic)
print("Valor p:", p_value)

t_statistic, p_value = stats.ttest_ind(random_minus_seq_sum, random_random_minus_sum)
print("seq + aleatorias")
print("Estadística de prueba (t):", t_statistic)
print("Valor p:", p_value)

t_statistic, p_value = stats.ttest_ind(random_plus_seq_sum,random_random_plus_sum)
print("seq - aleatorias")
print("Estadística de prueba (t):", t_statistic)
print("Valor p:", p_value)




def t_test_selection(sample1, sample2, alpha=0.05):
    # Prueba de igualdad de varianzas
    _, p_value = stats.levene(sample1, sample2)
    equal_var = p_value > alpha

    # Prueba de normalidad
    _, p_value = stats.shapiro(sample1)
    normal1 = p_value > alpha

    _, p_value = stats.shapiro(sample2)
    normal2 = p_value > alpha

    # Determinar qué prueba utilizar
    if normal1 and normal2:
        if equal_var:
            test_type = "Student's t-test"
        else:
            test_type = "Welch's t-test"
    else:
        test_type = "Non-parametric test (e.g., Mann-Whitney U test)"

    return test_type


# test_type = t_test_selection(plus_seq_sum, random_plus_seq_sum)
# print("Prueba recomendada:", test_type)






