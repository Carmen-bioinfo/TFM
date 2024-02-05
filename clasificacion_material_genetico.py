# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 17:55:56 2023

@author: Windows
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Cargamos nuestra tabla filtrada
input_path=r"C:\Users\Windows\OneDrive\Escritorio\TFM\MIS PROGRAMAS\SRX\tabla_filtrada_hsa.csv"
data_filtro=pd.read_table(input_path,sep='\t')
data_filtro=data_filtro.rename(columns={'id': 'Accession'})

#---------Cargamos la tabla donde los virus están anotados como host:human
input_path2=r"C:\Users\Windows\OneDrive\Escritorio\TFM\MIS PROGRAMAS\violinplot\clasificacion_por_material_genetico\sequences_human.csv"
data_NCBI=pd.read_table(input_path2, sep=',')

#----------fusionamos ambas tablas por las filas que tengan en común
fusionado = pd.merge(data_filtro, data_NCBI, on='Accession')

#----------pasos para ver qué filas faltan por añadir
fusion_filtro_human=data_filtro.merge(data_NCBI, on="Accession", how="left")
fusionado2_nulo = fusion_filtro_human[fusion_filtro_human['Molecule_type'].isnull()]
valores_accession = fusionado2_nulo['Accession'].tolist()
valores_separados = ', '.join(valores_accession)

#-----------añadimos las filas a la tabla de antes
input_path3=r"C:\Users\Windows\OneDrive\Escritorio\TFM\MIS PROGRAMAS\violinplot\clasificacion_por_material_genetico\sequences_null.csv"
data_NCBI_all=pd.read_table(input_path3, sep=',')

fusion_filtro_null=pd.merge(data_filtro, data_NCBI_all, on='Accession')
fusion_filtro_null = fusion_filtro_null.drop('Species', axis=1)
fusion_filtro_null = fusion_filtro_null.drop('USA', axis=1)

fusion_all=pd.concat([fusionado, fusion_filtro_null])
df_grouped = fusion_all.groupby(['Accession', 'Etiqueta', 'Molecule_type']).agg({'sum': 'sum'}).reset_index()

#-----------saber las filas que faltan por fusionar----------

merged_df = pd.merge(fusion_all, data_filtro, how='outer', indicator=True)
extra_rows = merged_df[merged_df['_merge'] == 'right_only']
extra_rows=extra_rows.drop('_merge', axis=1)
df_extra_rows= extra_rows.groupby(['Accession', 'Etiqueta']).agg({'sum': 'sum'}).reset_index()

#------------En estas filas hay que añadir manualmente el tipo de material genético, ya que no se encuentra en la base de datos NCBI virus---------

nuevos_valores = ['ssRNA(+)','ssDNA', 'ssRNA(+)', 'ssRNA(+)', 'ssRNA(+)', 'ssRNA(+)','ssRNA(+)','ssRNA(+)','ssRNA(+)','ssRNA(+)','ssRNA(+)','ssRNA(+)','ssRNA(+)','ssRNA(+)','ssRNA(+)','ssRNA(+)','ssRNA(+)','ssRNA(+)','ssRNA(+)','ssRNA(+)']
df_extra_rows['Molecule_type'] = nuevos_valores

#-------------lo fusionamos todo para obtener la tabla final---------------
fusion_final=pd.concat([df_grouped, df_extra_rows])
# fusion_final.to_string("fusion_final.txt")
#filas_nulas = fusion_final.loc[fusion_final['Etiqueta'].isnull()] #por si quisiéramos comprobar si vuelven a aparecer valores nulos, en este caso no


fusion_final['Molecule_type'].replace('ssDNA(+/-)', 'ssDNA', inplace=True)
fusion_final['Molecule_type'].replace('ssDNA(-)', 'ssDNA', inplace=True)
fusion_final['Molecule_type'].replace('unknown', 'circularDNA', inplace=True)

fusion_final["suma_log"] = np.log10(fusion_final["sum"]+1)

fusion_final.to_csv("fusion_final.csv", sep=",")
# #-------------ejemplo por si queremos seleccionar y ver los virus por tipo de material genético-----------
seleccion=fusion_final[fusion_final['Molecule_type']=='dsDNA']
seleccion = seleccion.sort_values(by='sum', ascending=False)

seleccion.to_string('dsDNA_abundancia.txt')

# #--------------pasamos a la parte del histograma---------------

# # Calcular la suma de la columna 'sum2' para cada etiqueta

suma_sum = fusion_final.groupby('Molecule_type')['sum'].sum()


dict_log={}

for indice, valor in suma_sum.iteritems():
    valor_log=np.log10(valor+1)
    dict_log[indice]=valor_log


dict_log=dict(sorted(dict_log.items(), key=lambda x: x[1], reverse=True))

# Crear el gráfico de barras con la suma de 'sum2' en el eje Y

# Obtener las claves y valores del diccionario
claves = dict_log.keys()
valores = dict_log.values()


# Crear el histograma
plt.bar(claves, valores)

# suma_sum.plot(kind='bar')


# Configurar el título y etiquetas de los ejes

plt.ylabel('log10(Abundancia total)')

plt.xticks(rotation=90)
plt.yticks(rotation=90)

# # Mostrar el gráfico de barras
plt.show()


#-----histogramas por desglose------
# virus_dsDNA=fusion_final[fusion_final["Molecule_type"] == "dsDNA"]

# virus_dsDNA=virus_dsDNA.sort_values(by="sum", ascending=False)

# virus_dsDNA=virus_dsDNA.head(15)

# virus_dsDNA['sum_log'] = np.log10(virus_dsDNA['sum'] + 1)
# # Eliminar ":human_virus_host" al final de los valores en la columna
# virus_dsDNA['Etiqueta'] = virus_dsDNA['Etiqueta'].str.replace(':human_virus_host', '')


# plt.figure(figsize=(10, 6))


# plt.bar(virus_dsDNA['Etiqueta'], virus_dsDNA['sum'])

# # Ajustar el espacio entre las etiquetas del eje X para evitar solapamientos
# plt.xticks(rotation=90)
# plt.ylabel('Abundancia total')
# plt.clean()


# #-----Segundo histograma por desglose------
"""
virus_ssRNA_RT=fusion_final[fusion_final["Molecule_type"] == "ssRNA-RT"]

virus_ssRNA_RT=virus_ssRNA_RT.sort_values(by="sum", ascending=False)

virus_ssRNA_RT=virus_ssRNA_RT.head(15)

virus_ssRNA_RT['sum_log'] = np.log10(virus_ssRNA_RT['sum'] + 1)
# Eliminar ":human_virus_host" al final de los valores en la columna
virus_ssRNA_RT['Etiqueta'] = virus_ssRNA_RT['Etiqueta'].str.replace(':human_virus_host', '')


plt.figure(figsize=(12, 6))


plt.bar(virus_ssRNA_RT['Etiqueta'], virus_ssRNA_RT['sum'])

# Ajustar el espacio entre las etiquetas del eje X para evitar solapamientos
plt.xticks(rotation=90)
plt.ylabel('Abundancia total')
plt.show()


"""