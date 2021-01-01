# Handling chromatine configuration data

import pandas as pd
import numpy as np


pd.options.mode.chained_assignment = None

# import data
df = pd.read_excel('4Dgenome.xlsx', header=0, parse_cols=[6, 7, 13])
df['Contact_Frequency'] = np.nan_to_num(df['Contact_Frequency'])
core_genes = pd.read_excel('Core Genes.xlsx', header=None)
non_core_genes = pd.read_excel('Non-core Genes.xlsx', header=None)


# data for core and non-core genes
core_data = []
for i in range(len(core_genes)):
    core_data.append(df['Agene'].str.find(core_genes.iloc[i, 0]))

non_core_data = []
for i in range(len(non_core_genes)):
    non_core_data.append(df['Agene'].str.find(non_core_genes.iloc[i, 0]))

# number of contacts for each core gene
sum_list_core = [[] for _ in range(65)]
for i in range(65):
    for j in range(len(core_data[1])):
        if core_data[i][j] >= 0:
            sum_list_core[i].append(j)
core_sums = [len(sum_list_core[i]) for i in range(len(sum_list_core))]

# number of contacts for each non core gene
sum_list_non_core = [[] for _ in range(len(non_core_genes))]
for i in range(len(non_core_genes)):
    for j in range(len(non_core_data[1])):
        if non_core_data[i][j] >= 0:
            sum_list_non_core[i].append(j)
non_core_sums = [len(sum_list_non_core[i]) for i in range(len(sum_list_non_core))]

# contact frequency for each core gene
core_contact = [[] for _ in range(65)]
for i in range(65):
    for j in range(len(core_data[1])):
        if core_data[i][j] is True:
            core_contact[i].append(df['Contact_Frequency'][j])
sum_core_contact = [sum(core_contact[i]) for i in range(len(core_contact))]

# contact frequency for each non-core gene
non_core_contact = [[] for _ in range(len(non_core_genes))]
for i in range(len(non_core_contact)):
    for j in range(len(non_core_data[1])):
        if non_core_data[i][j] is True:
            non_core_contact[i].append(df['Contact_Frequency'][j])
sum_non_core_contact = [sum(non_core_contact[i]) for i in range(len(non_core_contact))]

# Non-gene area contact point for core genes
core_nan = [[] for _ in range(65)]
for i in range(65):
    for j in range(len(core_data[1])):
        if core_data[i][j] is True:
            if pd.isnull(df['Bgene'][j]) is True:  # gives the non-gene regions
                core_nan[i].append(1)
sum_core_nan = [sum(core_nan[i]) for i in range(len(core_nan))]

# Non-gene area contact point for non-core genes
non_core_nan = [[] for _ in range(len(non_core_genes))]
for i in range(len(non_core_genes)):
    for j in range(len(non_core_data[1])):
        if non_core_data[i][j] is True:
            if pd.isnull(df['Bgene'][j]) is True:
                non_core_nan[i].append(1)
sum_non_core_nan = [sum(non_core_nan[i]) for i in range(len(non_core_nan))]

# search for core genes in 4D database
# core_4D = np.ndarray
# for ilist in range(len(core_data)):
#     for iseries in range(len(core_data[1])):
#         if core_data[ilist][iseries] == 0:
#             core_4D.append((df.loc[iseries].get_values()))


# writer = pd.ExcelWriter('Core 4D-genome.xlsx')
# core_4D.to_excel(writer, index_label=['Agene', 'Bg[ne', 'Contact_Frequency'], na_rep='NaN')

