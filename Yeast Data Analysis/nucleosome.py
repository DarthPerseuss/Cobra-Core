# Handling nucleosome modification data for core and non-core genes

import pandas as pd
import numpy as np
import statistics as st

df = pd.read_excel('YNA_Genome_Profiles.xlsx', header=0, parse_cols=[0, 3, 4, 5])
df_core = pd.read_excel('Core Genes.xlsx', header=None)
df_non_core = pd.read_excel('Non-core Genes.xlsx', header=None)


def match(target, string):
    """
    Function that evaluates if two strings are
    exactly equal.
    :param target: the target string
    :param string: the string to be compared to the target
    :return: bool
    """
    init = []
    if len(target) != len(string):
        return False
    for char in range(len(target)):
            if target[char] == string[char]:
                init.append(1)
            else:
                init.append(0)
    if sum(init) == len(target):
        return True
    else:
        return False

# find the matches and data corresponding to the nucleosome database for core genes
acore = []
for i in range(len(df)):
    for j in range(len(df_core)):
        if match(df.iloc[i, 0], (df_core.iloc[j, 0])) is True:
            acore.append(df.iloc[i].get_values())

# find the matches and data corresponding to the nucleosome database for non-core genes
anon_core = []
for i in range(len(df)):
    for j in range(len(df_non_core)):
        if match(df.iloc[i, 0], (df_non_core.iloc[j, 0])) is True:
            anon_core.append(df.iloc[i].get_values())

# find the matches and data corresponding to the nucleosome database for non-core genes


# **** After the previous procedure, it is time
# to import the extracted data and compare them
# with the histone modification data at hand; so
# this part is dependent on the previous part ,
# but also independent in its implementation****

df1 = pd.read_excel('Sequence_H3Pok.xlsx', skiprows=[i for i in range(31)], header=None)
df_core_pos = pd.read_excel('Core Gene Positions.xlsx', header=0)
df_non_core_pos = pd.read_excel('Non-core Gene Positions.xlsx', header=0)

# this block of code finds the core-related data from nucleosome modification data

histones = [[] for _ in range(len(df_core_pos))]
for i in range(len(df_core_pos)):
    for j in range(len(df1)):
        if min(df_core_pos['Coding Start'][i], df_core_pos['Coding End'][i]) <= df1[1][j] <= \
                max(df_core_pos['Coding Start'][i], df_core_pos['Coding End'][i]) and \
                df1[0][j] == df_core_pos['Chromosome'][i]:
            histones[i].append(df1.iloc[j].get_values())

# turn the acquired data into a dataframe
df_histones = pd.DataFrame(histones[0])
for i in range(1, len(histones)):
    df_histones = df_histones.append(pd.DataFrame(histones[i]))

df_histones.to_excel('Core Histone Mods.xlsx', index=False, header=['Chromosome', 'Position', 'H3K4me1 vs H3 (YPD)',
                                                                    'H3K4me2 vs H3 (YPD)', 'H3K4me3 vs H3 (YPD)',
                                                                    'H3K9ac vs H3 (YPD)', 'H3K14ac vs H3 (YPD)',
                                                                    'H3K14ac vs H3 (H2O2)', 'H3K14ac vs WCE (YPD)',
                                                                    'H3K36me3 vs H3 (YPD)',
                                                                    'H3K79me3 vs H3 (YPD)', 'H4 vs WCE (YPD)',
                                                                    'H4ac vs H3 (YPD)', 'H4ac vs H3 (H2O2)',
                                                                    'H3 vs WCE (YPD)', 'H3 vs WCE (H2O2)',
                                                                    'Esa1 vs WCE (YPD)', 'Gcn4 vs WCE (AA)',
                                                                    'Gcn5 vs WCE (YPD)', 'IgG control (YPD)',
                                                                    'No antibody control (YPD)'])

# same as above for non-core data
histones_nc = [[] for _ in range(len(df_non_core_pos))]
for i in range(len(df_non_core_pos)):
    for j in range(len(df1)):
        if min(df_non_core_pos['Coding Start'][i], df_non_core_pos['Coding End'][i]) <= df1[1][j] <= \
                max(df_non_core_pos['Coding Start'][i], df_non_core_pos['Coding End'][i]) and \
                df1[0][j] == df_non_core_pos['Chromosome'][i]:
            histones_nc[i].append(df1.iloc[j].get_values())

df_histones_nc = pd.DataFrame(histones_nc[0])
for i in range(1, len(histones_nc)):
    df_histones_nc = df_histones_nc.append(pd.DataFrame(histones_nc[i]))

df_histones_nc.to_excel('Non-core Histone Mods.xlsx', index=False, header=['Chromosome', 'Position',
                                                                           'H3K4me1 vs H3 (YPD)',
                                                                           'H3K4me2 vs H3 (YPD)', 'H3K4me3 vs H3 (YPD)',
                                                                           'H3K9ac vs H3 (YPD)', 'H3K14ac vs H3 (YPD)',
                                                                           'H3K14ac vs H3 (H2O2)',
                                                                           'H3K14ac vs WCE (YPD)',
                                                                           'H3K36me3 vs H3 (YPD)',
                                                                           'H3K79me3 vs H3 (YPD)', 'H4 vs WCE (YPD)',
                                                                           'H4ac vs H3 (YPD)', 'H4ac vs H3 (H2O2)',
                                                                           'H3 vs WCE (YPD)', 'H3 vs WCE (H2O2)',
                                                                           'Esa1 vs WCE (YPD)', 'Gcn4 vs WCE (AA)',
                                                                           'Gcn5 vs WCE (YPD)', 'IgG control (YPD)',
                                                                           'No antibody control (YPD)'])


# Acquiring the data for nucleosome ubiquitilation
df_ubiquitin = pd.read_excel('Sequence_H3Sch_Ubi.xlsx', skiprows=[i for i in range(26)], header=None)

histones_ubi = [[] for _ in range(len(df_core_pos))]
for i in range(len(df_core_pos)):
    for j in range(len(df_ubiquitin)):
        if min(df_core_pos['Coding Start'][i], df_core_pos['Coding End'][i]) <= df_ubiquitin[1][j] <= \
                max(df_core_pos['Coding Start'][i], df_core_pos['Coding End'][i]) and \
                df_ubiquitin[0][j] == df_core_pos['Chromosome'][i]:
            histones_ubi[i].append(df_ubiquitin.iloc[j].get_values())

df_histones_ubi = pd.DataFrame(histones_ubi[0])
for i in range(1, len(histones_ubi)):
    df_histones_ubi = df_histones_ubi.append(pd.DataFrame(histones_ubi[i]))

df_histones_ubi.to_excel('Core Histone Ubi.xlsx', index=False, header=['Chromosome', 'Position',
                                                                       'H2BK123ub vs mock', 'H3K4me3 vs input',
                                                                       'H3K36me3 vs input', 'H3K79me2 vs input',
                                                                       'H3K79me3 vs input', 'Rpb3 vs input',
                                                                       'H2BK123ub ubp8delta vs input',
                                                                       'H3K4me3 ubp8delta vs input',
                                                                       'H3K36me3 ubp8delta vs input',
                                                                       'H3K79me3 ubp8delta vs input',
                                                                       'Rpb3 ubp8delta vs input',
                                                                       ' H2BK123ub ubp10delta vs input',
                                                                       'H3K4me3 ubp10delta vs input',
                                                                       'H3K79me3 ubp10delta vs input'])

histones_ubi_nc = [[] for _ in range(len(df_non_core_pos))]
for i in range(len(df_non_core_pos)):
    for j in range(len(df_ubiquitin)):
        if min(df_non_core_pos['Coding Start'][i], df_non_core_pos['Coding End'][i]) <= df_ubiquitin[1][j] <= \
                max(df_non_core_pos['Coding Start'][i], df_non_core_pos['Coding End'][i]) and \
                df_ubiquitin[0][j] == df_non_core_pos['Chromosome'][i]:
            histones_ubi_nc[i].append(df_ubiquitin.iloc[j].get_values())


df_histones_ubi_nc = pd.DataFrame(histones_ubi_nc[0])
for i in range(1, len(histones_ubi_nc)):
    df_histones_ubi = df_histones_ubi.append(pd.DataFrame(histones_ubi_nc[i]))

df_histones_ubi_nc.to_excel('Non-core Histone Ubi.xlsx', index=False, header=['Chromosome', 'Position',
                                                                              'H2BK123ub vs mock', 'H3K4me3 vs input',
                                                                              'H3K36me3 vs input', 'H3K79me2 vs input',
                                                                              'H3K79me3 vs input', 'Rpb3 vs input',
                                                                              'H2BK123ub ubp8delta vs input',
                                                                              'H3K4me3 ubp8delta vs input',
                                                                              'H3K36me3 ubp8delta vs input',
                                                                              'H3K79me3 ubp8delta vs input',
                                                                              'Rpb3 ubp8delta vs input',
                                                                              'H2BK123ub ubp10delta vs input',
                                                                              'H3K4me3 ubp10delta vs input',
                                                                              'H3K79me3 ubp10delta vs input'])



# histones_ubi = [[] for _ in range(2)]
# for i in range(2):
#     for j in range(1, len(df_ubiquitin)):
#         flag = True
#         if flag is True:
#             if min(df_core_pos['Coding Start'][i], df_core_pos['Coding End'][i]) <= df_ubiquitin[1][j] <= \
#                             max(df_core_pos['Coding Start'][i], df_core_pos['Coding End'][i]) and \
#                             df_ubiquitin[0][j] == df_core_pos['Chromosome'][i]:
#                 histones_ubi[i].append(df_ubiquitin.iloc[j].get_values())
#             if df_ubiquitin[0][j - 1] != df_ubiquitin[0][j]:
#                 flag = False
#                 break

# This block of codes designates each position by its corresponding gene
df_core_pos = pd.read_excel('Core Gene Positions.xlsx', header=0)
core_his_ubi = pd.read_excel('Core Histone Ubi.xlsx', header=0)
i = 0
core_his_ubi.insert(16, 'Gene', 0)
for j in range(len(core_his_ubi)):
    if abs(core_his_ubi['Position'][j] - core_his_ubi['Position'][j+1]) < 100:
        core_his_ubi['Gene'][j] = df_core_pos['Systematic Name'][i]
    else:
        core_his_ubi['Gene'][j] = df_core_pos['Systematic Name'][i]
        core_his_ubi['Gene'][j+1] = df_core_pos['Systematic Name'][i+1]
        if i < 64:
            i += 1

# This block of code is added to address the mistake of choosing excel to
# read a large file. read_csv should be used if the process is to be re-implemented.
chrXVI = pd.read_csv('Sequence_H3Sch.csv', skiprows=[i for i in range(1048576)], header=None)
histones_xvi = [[] for _ in range(8)]
for i in range(59, 67):
    for j in range(len(chrXVI)):
        if min(df_core_pos['Coding Start'][i], df_core_pos['Coding End'][i]) <= chrXVI[1][j] <= \
                max(df_core_pos['Coding Start'][i], df_core_pos['Coding End'][i]) and \
                chrXVI[0][j] == df_core_pos['Chromosome'][i]:
            histones_xvi[i-59].append(chrXVI.iloc[j].get_values())

df_chrXVI = pd.DataFrame(histones_xvi[0])
for i in range(1, len(histones_xvi)):
    df_chrXVI = df_chrXVI.append(pd.DataFrame(histones_xvi[i]))

df_core_ubi = pd.read_excel('Core Histone Ubi.xlsx', header=0)
df_chrXVI.rename(inplace=True, columns={0: 'Chromosome', 1: 'Position', 2: 'H2BK123ub vs mock', 3: 'H3K4me3 vs input',
                                        4: 'H3K36me3 vs input', 5: 'H3K79me2 vs input', 6: 'H3K79me3 vs input',
                                        7: 'Rpb3 vs input', 8: 'H2BK123ub ubp8delta vs input',
                                        9: 'H3K4me3 ubp8delta vs input', 10: 'H3K36me3 ubp8delta vs input',
                                        11: 'H3K79me3 ubp8delta vs input', 12: 'Rpb3 ubp8delta vs input',
                                        13: 'H2BK123ub ubp10delta vs input', 14: 'H3K4me3 ubp10delta vs input',
                                        15: 'H3K79me3 ubp10delta vs input'})

# Does not result in a clean dataframe, maybe look at later
a = pd.concat([df_core_ubi, df_chrXVI], ignore_index=True, join='inner', verify_integrity=True)

# Manually copy-paste to the original excel file
df_chrXVI.to_excel('#NAME.xlsx', index=False, header=['Chromosome', 'Position',
                                                      'H2BK123ub vs mock', 'H3K4me3 vs input',
                                                      'H3K36me3 vs input', 'H3K79me2 vs input',
                                                      'H3K79me3 vs input', 'Rpb3 vs input',
                                                      'H2BK123ub ubp8delta vs input',
                                                      'H3K4me3 ubp8delta vs input',
                                                      'H3K36me3 ubp8delta vs input',
                                                      'H3K79me3 ubp8delta vs input',
                                                      'Rpb3 ubp8delta vs input',
                                                      'H2BK123ub ubp10delta vs input',
                                                      'H3K4me3 ubp10delta vs input',
                                                      'H3K79me3 ubp10delta vs input'])

# assigning genes for non-core ubiquitilation
df_non_core_pos = pd.read_excel('Non-core Gene Positions.xlsx', header=0)
ncore_his_ubi = pd.read_excel('Non-core Histone Ubi.xlsx', header=0)
i = 0
ncore_his_ubi.insert(16, 'Gene', 0)
for j in range(len(ncore_his_ubi)):
    if abs(ncore_his_ubi['Position'][j] - ncore_his_ubi['Position'][j+1]) < 100:
        ncore_his_ubi['Gene'][j] = df_non_core_pos['Systematic Name'][i]
    else:
        ncore_his_ubi['Gene'][j] = df_non_core_pos['Systematic Name'][i]
        ncore_his_ubi['Gene'][j+1] = df_non_core_pos['Systematic Name'][i+1]
        if i < len(df_non_core_pos) - 1:
            i += 1
# same as above for non-core genes
chrXVI = pd.read_csv('Sequence_H3Sch.csv', skiprows=[i for i in range(1048576)], header=None)
histones_xvi = [[] for _ in range(71)]
for i in range(761, 832):
    for j in range(len(chrXVI)):
        if min(df_non_core_pos['Coding Start'][i], df_non_core_pos['Coding End'][i]) <= chrXVI[1][j] <= \
                max(df_non_core_pos['Coding Start'][i], df_non_core_pos['Coding End'][i]) and \
                chrXVI[0][j] == df_non_core_pos['Chromosome'][i]:
            histones_xvi[i-761].append(chrXVI.iloc[j].get_values())

df_chrXVI_nc = pd.DataFrame(histones_xvi[0])
for i in range(1, len(histones_xvi)):
    df_chrXVI_nc = df_chrXVI_nc.append(pd.DataFrame(histones_xvi[i]))

df_chrXVI_nc.rename(inplace=True,
                    columns={0: 'Chromosome', 1: 'Position', 2: 'H2BK123ub vs mock', 3: 'H3K4me3 vs input',
                             4: 'H3K36me3 vs input', 5: 'H3K79me2 vs input', 6: 'H3K79me3 vs input',
                             7: 'Rpb3 vs input', 8: 'H2BK123ub ubp8delta vs input',
                             9: 'H3K4me3 ubp8delta vs input', 10: 'H3K36me3 ubp8delta vs input',
                             11: 'H3K79me3 ubp8delta vs input', 12: 'Rpb3 ubp8delta vs input',
                             13: 'H2BK123ub ubp10delta vs input', 14: 'H3K4me3 ubp10delta vs input',
                             15: 'H3K79me3 ubp10delta vs input'})

df_chrXVI_nc.to_excel('#NAME.xlsx', index=False, header=['Chromosome', 'Position',
                                                         'H2BK123ub vs mock', 'H3K4me3 vs input',
                                                         'H3K36me3 vs input', 'H3K79me2 vs input',
                                                         'H3K79me3 vs input', 'Rpb3 vs input',
                                                         'H2BK123ub ubp8delta vs input',
                                                         'H3K4me3 ubp8delta vs input',
                                                         'H3K36me3 ubp8delta vs input',
                                                         'H3K79me3 ubp8delta vs input',
                                                         'Rpb3 ubp8delta vs input',
                                                         'H2BK123ub ubp10delta vs input',
                                                         'H3K4me3 ubp10delta vs input',
                                                         'H3K79me3 ubp10delta vs input'])

# assigning genes to core histone modifications
df_core_pos = pd.read_excel('Core Gene Positions.xlsx', header=0)
core_his_mod = pd.read_excel('Core Histone Mods.xlsx', header=0)
core_his_mod.insert(21, 'Gene', 0)
for i in range(len(df_core_pos)):
    for j in range(len(core_his_mod)):
        if min(df_core_pos['Coding Start'][i], df_core_pos['Coding End'][i]) <= core_his_mod['Position'][j] <= \
            max(df_core_pos['Coding Start'][i], df_core_pos['Coding End'][i]) and \
                core_his_mod['Chromosome'][j] == df_core_pos['Chromosome'][i]:
                core_his_mod['Gene'][j] = df_core_pos['Systematic Name'][i]

core_his_mod.to_excel('Core Histone Mods.xlsx', index=False, header=['Chromosome', 'Position', 'H3K4me1 vs H3 (YPD)',
                                                                     'H3K4me2 vs H3 (YPD)', 'H3K4me3 vs H3 (YPD)',
                                                                     'H3K9ac vs H3 (YPD)', 'H3K14ac vs H3 (YPD)',
                                                                     'H3K14ac vs H3 (H2O2)', 'H3K14ac vs WCE (YPD)',
                                                                     'H3K36me3 vs H3 (YPD)',
                                                                     'H3K79me3 vs H3 (YPD)', 'H4 vs WCE (YPD)',
                                                                     'H4ac vs H3 (YPD)', 'H4ac vs H3 (H2O2)',
                                                                     'H3 vs WCE (YPD)', 'H3 vs WCE (H2O2)',
                                                                     'Esa1 vs WCE (YPD)', 'Gcn4 vs WCE (AA)',
                                                                     'Gcn5 vs WCE (YPD)', 'IgG control (YPD)',
                                                                     'No antibody control (YPD)', 'Gene'])

# assigning genes to non-core histone modifications
df_non_core_pos = pd.read_excel('Non-core Gene Positions.xlsx', header=0)
ncore_his_mod = pd.read_excel('Non-core Histone Mods.xlsx', header=0)
ncore_his_mod.insert(21, 'Gene', 0)
for i in range(len(df_non_core_pos)):
    for j in range(len(ncore_his_mod)):
        if min(df_non_core_pos['Coding Start'][i], df_non_core_pos['Coding End'][i]) <= ncore_his_mod['Position'][j] <= \
            max(df_non_core_pos['Coding Start'][i], df_non_core_pos['Coding End'][i]) and \
                ncore_his_mod['Chromosome'][j] == df_non_core_pos['Chromosome'][i]:
                ncore_his_mod['Gene'][j] = df_non_core_pos['Systematic Name'][i]

ncore_his_mod.to_excel('Non-core Histone Mods.xlsx', index=False,
                       header=['Chromosome', 'Position', 'H3K4me1 vs H3 (YPD)',
                               'H3K4me2 vs H3 (YPD)', 'H3K4me3 vs H3 (YPD)',
                               'H3K9ac vs H3 (YPD)', 'H3K14ac vs H3 (YPD)',
                               'H3K14ac vs H3 (H2O2)', 'H3K14ac vs WCE (YPD)',
                               'H3K36me3 vs H3 (YPD)',
                               'H3K79me3 vs H3 (YPD)', 'H4 vs WCE (YPD)',
                               'H4ac vs H3 (YPD)', 'H4ac vs H3 (H2O2)',
                               'H3 vs WCE (YPD)', 'H3 vs WCE (H2O2)',
                               'Esa1 vs WCE (YPD)', 'Gcn4 vs WCE (AA)',
                               'Gcn5 vs WCE (YPD)', 'IgG control (YPD)',
                               'No antibody control (YPD)', 'Gene'])


# This block of code separates entries based on the corresponding gene
# and computes means for each variable in the original database.
core_his_mod = pd.read_excel('Core Histone Mods.xlsx', header=0)
df_core_pos = pd.read_excel('Core Gene Positions.xlsx', header=0)
a = [[] for _ in range(len(df_core_pos))]
for i in range(len(a)):
    a[i].append(core_his_mod[core_his_mod['Gene'] == df_core_pos['Systematic Name'][i]])

b = [[] for _ in range(len(a))]
for i in range(len(a)):
    for j in range(2,21):
        b[i].append(st.mean(a[i][0].iloc[:, j]))


for i in range(len(b)):
    b[i] = np.asarray(b[i])
df_histones_mod_mean = pd.DataFrame(columns={'H3K4me1 vs H3 (YPD)',
                                             'H3K4me2 vs H3 (YPD)', 'H3K4me3 vs H3 (YPD)',
                                             'H3K9ac vs H3 (YPD)', 'H3K14ac vs H3 (YPD)',
                                             'H3K14ac vs H3 (H2O2)', 'H3K14ac vs WCE (YPD)',
                                             'H3K36me3 vs H3 (YPD)',
                                             'H3K79me3 vs H3 (YPD)', 'H4 vs WCE (YPD)',
                                             'H4ac vs H3 (YPD)', 'H4ac vs H3 (H2O2)',
                                             'H3 vs WCE (YPD)', 'H3 vs WCE (H2O2)',
                                             'Esa1 vs WCE (YPD)', 'Gcn4 vs WCE (AA)',
                                             'Gcn5 vs WCE (YPD)', 'IgG control (YPD)',
                                             'No antibody control (YPD)'}, index=[i for i in range(len(a))])
for i in range(len(b)):
    for j in range(len(b[0])):
        df_histones_mod_mean.iloc[i, j] = b[i][j]

df_histones_mod_mean.to_excel('Core Histone Mods.xlsx', index=False, sheet_name='sheet2',
                              header=['H3K4me1 vs H3 (YPD)',
                                      'H3K4me2 vs H3 (YPD)', 'H3K4me3 vs H3 (YPD)',
                                      'H3K9ac vs H3 (YPD)', 'H3K14ac vs H3 (YPD)',
                                      'H3K14ac vs H3 (H2O2)', 'H3K14ac vs WCE (YPD)',
                                      'H3K36me3 vs H3 (YPD)',
                                      'H3K79me3 vs H3 (YPD)', 'H4 vs WCE (YPD)',
                                      'H4ac vs H3 (YPD)', 'H4ac vs H3 (H2O2)',
                                      'H3 vs WCE (YPD)', 'H3 vs WCE (H2O2)',
                                      'Esa1 vs WCE (YPD)', 'Gcn4 vs WCE (AA)',
                                      'Gcn5 vs WCE (YPD)', 'IgG control (YPD)',
                                      'No antibody control (YPD)'])

# Same as above for non-core genes
ncore_his_mod = pd.read_excel('Non-core Histone Mods.xlsx', header=0)
df_non_core_pos = pd.read_excel('Non-core Gene Positions.xlsx', header=0)
a = [[] for _ in range(len(df_non_core_pos))]
for i in range(len(a)):
    a[i].append(ncore_his_mod[ncore_his_mod['Gene'] == df_non_core_pos['Systematic Name'][i]])

b = [[] for _ in range(len(a))]
for i in range(len(a)):
    for j in range(2, 21):
        try:  # this is to avoid errors and ensure the continuation of implementation
            b[i].append(st.mean(a[i][0].iloc[:, j]))
        except st.StatisticsError:
            b[i].append('NaN')

for i in range(len(b)):
    b[i] = np.asarray(b[i])
df_histones_mod_mean_nc = pd.DataFrame(columns={'H3K4me1 vs H3 (YPD)',
                                                'H3K4me2 vs H3 (YPD)', 'H3K4me3 vs H3 (YPD)',
                                                'H3K9ac vs H3 (YPD)', 'H3K14ac vs H3 (YPD)',
                                                'H3K14ac vs H3 (H2O2)', 'H3K14ac vs WCE (YPD)',
                                                'H3K36me3 vs H3 (YPD)',
                                                'H3K79me3 vs H3 (YPD)', 'H4 vs WCE (YPD)',
                                                'H4ac vs H3 (YPD)', 'H4ac vs H3 (H2O2)',
                                                'H3 vs WCE (YPD)', 'H3 vs WCE (H2O2)',
                                                'Esa1 vs WCE (YPD)', 'Gcn4 vs WCE (AA)',
                                                'Gcn5 vs WCE (YPD)', 'IgG control (YPD)',
                                                'No antibody control (YPD)'}, index=[i for i in range(len(a))])
for i in range(len(b)):
    for j in range(len(b[0])):
        df_histones_mod_mean_nc.iloc[i, j] = b[i][j]

df_histones_mod_mean_nc = df_histones_mod_mean_nc[df_histones_mod_mean_nc['H3K4me1 vs H3 (YPD)'] != 'NaN']

writer = pd.ExcelWriter('Non-core Histone Mods Mean.xlsx')
df_histones_mod_mean_nc.to_excel(writer, index=False, sheet_name='sheet2',
                                 header=['H3K4me1 vs H3 (YPD)',
                                         'H3K4me2 vs H3 (YPD)', 'H3K4me3 vs H3 (YPD)',
                                         'H3K9ac vs H3 (YPD)', 'H3K14ac vs H3 (YPD)',
                                         'H3K14ac vs H3 (H2O2)', 'H3K14ac vs WCE (YPD)',
                                         'H3K36me3 vs H3 (YPD)',
                                         'H3K79me3 vs H3 (YPD)', 'H4 vs WCE (YPD)',
                                         'H4ac vs H3 (YPD)', 'H4ac vs H3 (H2O2)',
                                         'H3 vs WCE (YPD)', 'H3 vs WCE (H2O2)',
                                         'Esa1 vs WCE (YPD)', 'Gcn4 vs WCE (AA)',
                                         'Gcn5 vs WCE (YPD)', 'IgG control (YPD)',
                                         'No antibody control (YPD)'])
writer.save()

#
#
#
# ****************** Core Histone Ubi**********************

core_his_ubi = pd.read_excel('Core Histone Ubi.xlsx', header=0)
df_core_pos = pd.read_excel('Core Gene Positions.xlsx', header=0)
a = [[] for _ in range(len(df_core_pos))]
for i in range(len(a)):
    a[i].append(core_his_ubi[core_his_ubi['Gene'] == df_core_pos['Systematic Name'][i]])

b = [[] for _ in range(len(a))]
for i in range(len(a)):
    for j in range(2, 16):
        try:  # this is to avoid errors and ensure the continuation of implementation
            b[i].append(st.mean(a[i][0].iloc[:, j]))
        except st.StatisticsError:
            b[i].append('NaN')


for i in range(len(b)):
    b[i] = np.asarray(b[i])
df_histones_ubi_mean = pd.DataFrame(columns={'H2BK123ub vs mock', 'H3K4me3 vs input',
                                             'H3K36me3 vs input', 'H3K79me2 vs input',
                                             'H3K79me3 vs input', 'Rpb3 vs input',
                                             'H2BK123ub ubp8delta vs input',
                                             'H3K4me3 ubp8delta vs input',
                                             'H3K36me3 ubp8delta vs input',
                                             'H3K79me3 ubp8delta vs input',
                                             'Rpb3 ubp8delta vs input',
                                             'H2BK123ub ubp10delta vs input',
                                             'H3K4me3 ubp10delta vs input',
                                             'H3K79me3 ubp10delta vs input'}, index=[i for i in range(len(a))])
for i in range(len(b)):
    for j in range(len(b[0])):
        df_histones_ubi_mean.iloc[i, j] = b[i][j]

writer = pd.ExcelWriter('Core Histone Ubi Mean.xlsx')
df_histones_ubi_mean.to_excel(writer, index=False, sheet_name='sheet2',
                              header=('H2BK123ub vs mock', 'H3K4me3 vs input',
                                      'H3K36me3 vs input', 'H3K79me2 vs input',
                                      'H3K79me3 vs input', 'Rpb3 vs input',
                                      'H2BK123ub ubp8delta vs input',
                                      'H3K4me3 ubp8delta vs input',
                                      'H3K36me3 ubp8delta vs input',
                                      'H3K79me3 ubp8delta vs input',
                                      'Rpb3 ubp8delta vs input',
                                      'H2BK123ub ubp10delta vs input',
                                      'H3K4me3 ubp10delta vs input',
                                      'H3K79me3 ubp10delta vs input'))
writer.save()

#
#
#
# ********* Non-core ubi mean ***********

ncore_his_ubi = pd.read_excel('Non-core Histone Ubi.xlsx', header=0)
df_ncore_pos = pd.read_excel('Non-core Gene Positions.xlsx', header=0)
a = [[] for _ in range(len(df_ncore_pos))]
for i in range(len(a)):
    a[i].append(ncore_his_ubi[ncore_his_ubi['Gene'] == df_ncore_pos['Systematic Name'][i]])

b = [[] for _ in range(len(a))]
for i in range(len(a)):
    for j in range(2, 16):
        try:  # this is to avoid errors and ensure the continuation of implementation
            b[i].append(st.mean(a[i][0].iloc[:, j]))
        except st.StatisticsError:
            b[i].append('NaN')

for i in range(len(b)):
    b[i] = np.asarray(b[i])
df_histones_ubi_nc_mean = pd.DataFrame(columns={'H2BK123ub vs mock', 'H3K4me3 vs input',
                                                'H3K36me3 vs input', 'H3K79me2 vs input',
                                                'H3K79me3 vs input', 'Rpb3 vs input',
                                                'H2BK123ub ubp8delta vs input',
                                                'H3K4me3 ubp8delta vs input',
                                                'H3K36me3 ubp8delta vs input',
                                                'H3K79me3 ubp8delta vs input',
                                                'Rpb3 ubp8delta vs input',
                                                'H2BK123ub ubp10delta vs input',
                                                'H3K4me3 ubp10delta vs input',
                                                'H3K79me3 ubp10delta vs input'}, index=[i for i in range(len(a))])
for i in range(len(b)):
    for j in range(len(b[0])):
        df_histones_ubi_nc_mean.iloc[i, j] = b[i][j]

writer = pd.ExcelWriter('Non-core Histone Ubi Mean.xlsx')
df_histones_ubi_nc_mean.to_excel(writer, index=False, sheet_name='sheet2',
                                 header=('H2BK123ub vs mock', 'H3K4me3 vs input',
                                         'H3K36me3 vs input', 'H3K79me2 vs input',
                                         'H3K79me3 vs input', 'Rpb3 vs input',
                                         'H2BK123ub ubp8delta vs input',
                                         'H3K4me3 ubp8delta vs input',
                                         'H3K36me3 ubp8delta vs input',
                                         'H3K79me3 ubp8delta vs input',
                                         'Rpb3 ubp8delta vs input',
                                         'H2BK123ub ubp10delta vs input',
                                         'H3K4me3 ubp10delta vs input',
                                         'H3K79me3 ubp10delta vs input'))
writer.save()
