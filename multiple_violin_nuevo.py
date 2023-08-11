# -*- coding: utf-8 -*-
"""
Created on Thu May 11 17:25:54 2023

@author: Windows
"""

import plotly.io as io
io.renderers.default='browser'
import plotly.graph_objects as go
import plotly.express as px


import numpy as np
import pandas as pd

# import argparse
# import Annotation

# parser = argparse.ArgumentParser()
# parser.add_argument("-f", "--file", help="matriz file")
# args = parser.parse_args()
# Annotation.kunitz_annotate_ndomains(args.file)

input_path=r"C:\Users\Windows\OneDrive\Escritorio\TFM\MIS PROGRAMAS\SRX\tabla_sin_filtrar_hsa.csv"

# El de abajo es C.elegans-Drosophila
#input_path=r"C:\Users\Windows\OneDrive\Escritorio\TFM\MIS PROGRAMAS\SRX\RESULTADOS\CONTROL NEGATIVO ELEGANS-DROSOPHILA\matriz_estadistica_etiqueta.csv"

#El de abajo es plantas
#input_path=r"C:\Users\Windows\OneDrive\Escritorio\TFM\MIS PROGRAMAS\SRX\RESULTADOS\CONTROL_NEGATIVO_PLANTAS\tabla_sin_filtrar_arab_maiz.csv"



# ESTA PARTE DEL CÓDIGO SOLO LEE Y PREPARA LOS DATOS PARA PODER HACER LA REPRESENTACIÓN
tabla = pd.read_csv(input_path, header=0, delimiter="\t")


# seleccion=tabla.head(20)

seleccion=tabla.copy()

split_genlib=seleccion['id'].str.split('_', n=2, expand=True)
nombre_NC=split_genlib[0]+"_"+split_genlib[1]
 
tipo_seq=split_genlib[2]

tipo_seq = tipo_seq.replace(to_replace=np.nan, value='seq')

del(seleccion["id"])

seleccion.loc[:, "tipo_seq"] = tipo_seq
seleccion.loc[:, "gene"] = nombre_NC

# Obtener el índice actual de las columnas
cols = seleccion.columns.tolist()

# Mover la columna "col3" al primer lugar
cols.insert(0, cols.pop(cols.index('gene')))
cols.insert(1, cols.pop(cols.index('tipo_seq')))

# # Reindexar el dataframe con el nuevo orden de columnas
df_trabajo = seleccion.reindex(columns=cols)

# Create a list of colors for each category
category_colors = {
    1: "blue",
    2: "red",
    3: "green",
    4: "purple"
}

data_violin = pd.DataFrame(dict(Score = df_trabajo["sum"], Group = df_trabajo["tipo_seq"]))

order = ["seq", "seq_aleatoria", "seq_reverso_complementaria", "seq_aleatoria_complementaria"]
data_violin["Group"] = pd.Categorical(data_violin["Group"], categories=order)

data_violin = data_violin.sort_values("Group")

data_violin["Suma_log"] = np.log10(data_violin["Score"]+1)

data_violin['GroupMod'] = pd.factorize(data_violin['Group'])[0] + 1



fig = px.violin(data_violin, y="Suma_log", x="GroupMod", color='GroupMod',
                box=True, points=False,
                hover_data=data_violin.columns)

rug_colors = [category_colors[category] for category in data_violin['GroupMod']]

# for valor in rug_colors:
#     if valor == 'red':
#         rug_colors[rug_colors.index(valor)] = 'blue'
#     elif valor == 'purple':
#         rug_colors[rug_colors.index(valor)] = 'green'



rug_trace = go.Scatter(x=data_violin['GroupMod']-0.3, y=data_violin['Suma_log'],
                        mode='markers',
                        marker=dict(color=rug_colors, size=10, symbol='line-ew-open'),
                        hoverinfo='skip', showlegend=False)

fig.add_trace(rug_trace)

fig.update_layout(
    xaxis=dict(
        tickmode='array',
        tickvals=[0.7,1,1.7,2,2.7,3,3.7,4],
        ticktext=["","Sequence plus","","Sequence plus random","","Sequence minus","","Sequence minus random"],
        tickangle=0
    )
)

fig.update_layout(showlegend=False)


fig.update_layout(yaxis_title="log10(suma de RPMs)")
fig.update_xaxes(title='')

# fig.update_layout(title={'text': "<b>CUMULATIVE ABUNDANCE OF VIRAL RNA AGAINST THEIR NEGATIVE CONTROLS</b>",
#         'x': 0.5,
#         'y': 0.95,
#         'xanchor': 'center',
#         'yanchor': 'top',
#         })

fig.show()