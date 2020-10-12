from __future__ import print_function
import cobra
import libsbml
from itertools import chain
import numpy as np
import gc

gc.enable()
# gc_thresh = gc.get_threshold()
# print(gc_thresh)


model = cobra.io.read_sbml_model('S.pombe genome model.xml')
# prints the objective reaction
print(model.reactions[len(model.reactions)-1])
print(model.objective.expression)
# sets the objective function
model.objective = model.reactions[len(model.reactions)-1]
# runs FBA
solution = model.optimize()

# to find reactions that have glucose in them
glucose = []
for i in range(len(model.reactions)):
    if 'GLC' in str(model.reactions[i]):
        glucose.append(i)

for react in glucose:
    print(model.reactions[react])

# doubly open reactions
rev_react = []
for i in range(len(model.reactions)):
    if model.reactions[i].bounds[0] != 0:
        rev_react.append(i)

react_comp = []
for i in range(len(model.reactions)):
    react_comp.append(model.reactions[i].get_compartments())
# un-nests the nested list
react_comp = list(chain.from_iterable(react_comp))
# outputs the unique members of the list
print(set(react_comp))

exc_rxns = model.exchanges
# uptake reactions are exchange reactions with negative lower bounds
removed = []
for i in range(len(exc_rxns)):
    if exc_rxns[i].bounds[0] >= 0:
        removed.append(exc_rxns[i])

upt_rxns = list(set(exc_rxns) - set(removed))
blocked_rxns = cobra.flux_analysis.find_blocked_reactions(model)

# in order to get more fidel results, we omit the blocked reactions
# from the set of initial uptake reactions
upt_rxnIDs = []
for i in range(len(upt_rxns)):
    upt_rxnIDs.append(upt_rxns[i]._id)

# this code gives the reactions as string IDs
finalUpt_rxns = list(set(upt_rxnIDs) - (set(blocked_rxns) & set(upt_rxnIDs)))

# to obtain the reaction numbers for further processing
fuRxn_nums = []
rxn_dict = {model.reactions[i]: i for i in range(len(model.reactions))}
for i in rxn_dict:
    for j in range(len(finalUpt_rxns)):
        if finalUpt_rxns[j] in str(i):
            fuRxn_nums.append((rxn_dict[i]))
        else:
            pass


exp_conds = 50000
solutions = [[] for _ in range(exp_conds)]
sample_models = [[] for _ in range(exp_conds)]
rand_LBs = np.random.uniform(-20, 0, (exp_conds, len(fuRxn_nums)))
active_rxns = np.zeros(shape=(exp_conds, len(model.reactions)))
copy = model.copy()


for row in range(exp_conds):
    sample_models[0].append(copy)
    sample_models[0][0].objective = model.reactions[len(model.reactions)-1]
    for col in range(len(fuRxn_nums)):
        sample_models[0][0].reactions[fuRxn_nums[col]].bounds = (rand_LBs[row][col], 1000)
    solutions[0].append(sample_models[0][0].optimize())
    # Here each model is deleted after the required computation is handled.
    # This is to free up memory using garbage collection. The list is popped
    # from the right-end side as the loop cycles.
    if row < exp_conds:
        del sample_models[0]
    for j in range(len(model.reactions)):
        if solutions[0][0].fluxes[j] > 0:
            active_rxns[row][j] = 1
    del solutions[0]
    gc.collect()

sum_activeRXNs = sum(active_rxns)
active_rxnNUMs = []
for i in range(len(sum_activeRXNs)):
    if sum_activeRXNs[i] == exp_conds:
        active_rxnNUMs.append(i)

init_core = [138,
             152,
             153,
             156,
             176,
             177,
             185,
             186,
             187,
             188,
             207,
             208,
             209,
             213,
             215,
             232,
             250,
             288,
             290,
             326,
             329,
             342,
             349,
             382,
             446,
             448,
             484,
             505,
             527,
             536,
             537,
             540,
             559,
             590,
             600,
             601,
             612,
             614,
             617,
             619,
             620,
             626,
             654,
             712,
             714,
             715,
             718,
             721,
             723,
             724,
             726,
             752,
             756,
             758,
             761,
             762,
             769,
             778,
             788,
             795,
             805,
             814,
             817,
             819,
             833,
             836,
             846,
             862,
             874,
             877,
             878,
             888,
             893,
             910,
             928,
             931,
             939,
             940,
             962,
             992,
             1033,
             1035,
             1050,
             1065,
             1069,
             1080,
             1102,
             1116,
             1194,
             1204,
             1240,
             1242,
             1285,
             1296,
             1302,
             1306,
             1312,
             1314,
             1316,
             1320,
             1339,
             1340,
             1353,
             1370,
             1371,
             1373,
             1374,
             1376,
             1384,
             1387,
             1406,
             1408,
             1410,
             1433,
             1445,
             1448,
             1449,
             1457,
             1480,
             1481,
             1486,
             1510,
             1517,
             1518,
             1537,
             1538,
             1539,
             1541,
             1548,
             1549,
             1583,
             1585,
             1596,
             1600,
             1633,
             1640,
             1645,
             1647,
             1685,
             1686,
             1687,
             1688,
             1689,
             1690,
             1691,
             1692]

# cobra.io.save_matlab_model(model, 'S_pombe.xml')
# CO2 : 516 ; ergosterol: 658 ; H2O : 922 ; ammonium: 1203 ;
# oxygen: 1240 ; phosphate: 1343 ; sulphate: 1523 ; zymosterol : 1693

# to find extracelluar transport reactions
for i in range(len(model.reactions)):
    if len(model.reactions[i].reactants) == 0:
        print("%s : %i" % (model.reactions[i], i))

mant_e = cobra.Metabolite('mant_e', name='Manninotriose', compartment='e')
gal_e = cobra.Metabolite('gal_e', name='Galactose', compartment='e')
H2O_e = cobra.Metabolite('H2O_e', name='Water', compartment='e')
melib_e = cobra.Metabolite('melib_e', name='Melibiose', compartment='e')
raffin_e = cobra.Metabolite('raffin_e', name='Raffinose', compartment='e')
sta_e = cobra.Metabolite('sta_e', name='Stachyose', compartment='e')
GALI_e = cobra.Metabolite('GALI_e', name='Galactinol', compartment='e')
inost_e = cobra.Metabolite('inost_e', name='myo-Inositol', compartment='e')
emp_e = cobra.Metabolite('emp_e', name='Epimelibioise', compartment='e')
man_e = cobra.Metabolite('man_e', name='D-Mannose', compartment='e')
ggl_e = cobra.Metabolite('ggl_e', name='Galactosylglycerol', compartment='e')
glyc_e = cobra.Metabolite('glyc_e', name='Glycerol', compartment='e')
melt_e = cobra.Metabolite('melt_e', name='Melibiitol', compartment='e')
sbtD_e = cobra.Metabolite('sbtD_e', name='D-sorbitol', compartment='e')
GDP_e = cobra.Metabolite('GDP_e', name='GDP', compartment='e')
GMP_e = cobra.Metabolite('GMP_e', name='GMP', compartment='e')
H_e = cobra.Metabolite('H_e', name='Hydrogen', compartment='e')
Pi_e = cobra.Metabolite('Pi_e', name='Phosphate', compartment='e')
dgala_e = cobra.Metabolite('dgala_e', name='Digalactosylceramide', compartment='e')
gala_e = cobra.Metabolite('gala_e', name='Galactosylceramide', compartment='e')
a_13BDglcn_e = cobra.Metabolite('a_13BDglcn_e', name='1,3-beta-D-Glucan', compartment='e')
GLC_e = cobra.Metabolite('GLC_e', name='D-Glucose', compartment='e')
a_2hb_e = cobra.Metabolite('a_2hb_e', name='2-Hydroxybutyrate', compartment='e')
zymst_e = cobra.Metabolite('zymst_e', name='Zymosterol', compartment='e')

# Adding the reaction for zymosterol reversible transport
zrt = cobra.Reaction('RXNZymst')
zrt.name = 'Zymosterol Reversible Transport'
zrt.lower_bound = -1000.0
zrt.upper_bound = 1000.0
# First the model has to be added and then the metabolites.
model.add_reaction(zrt)
zrt.add_metabolites({"zymst_c": 1}, reversibly=True)

# model.add_metabolites([mant_e, melib_e, H2O_e, gal_e])
# model.reactions[30].add_metabolites({"gal_e": 1.0, "melib_e": 1.0})

# Outputs transport reactions
revTransRXNs = []
for i in range(len(model.reactions)):
    if not model.reactions[i].reactants and model.reactions[i].compartments == {'cell'}:
        revTransRXNs.append(i)
        print(model.reactions[i], i)

# TO BE OMITTED FROM THE LIST ABOVE: fe2_c : 761, K_c : 1046;
[revTransRXNs.remove(i) for i in [516, 658, 922, 1203, 1240, 1343, 1523, 1693, 761, 1046]]
check_list = []

[check_list.append(1) if i in revTransRXNs else check_list.append(0)
 for i in [516, 658, 922, 1203, 1240, 1343, 1523, 1693, 761, 1046]]  # To make sure.


# To assign unlimited flux to the certain compounds required
# for a minimal uptake medium.
min_req = [516, 658, 922, 1203, 1240, 1343, 1523, 1693]
oneC_min_models = [[] for _ in range(len(revTransRXNs))]
for i in range(len(revTransRXNs)):
    oneC_min_models[i].append(model.copy())
    for j in min_req:
        oneC_min_models[i][0].reactions[j].lower_bound = -1000

# To assign fluxes to single-carbon source minimal uptake models
for i in range(len(oneC_min_models)):
    oneC_min_models[i][0].reactions[revTransRXNs[i]].lower_bound = -900

# To 'minimize' the models to being single carbon source.
for i in range(len(oneC_min_models)):
    for j in range(len(model.reactions)):
        if oneC_min_models[i][0].reactions[j].lower_bound == -900 or j in min_req:
            pass
        else:
            oneC_min_models[i][0].reactions[j].lower_bound = 0

min_activeRXNs = np.zeros(shape=(len(oneC_min_models), len(model.reactions)))
min_model_sols = [[] for _ in range(len(oneC_min_models))]

# to assign objective reactions to minimal models
for i in range(len(oneC_min_models)):
    oneC_min_models[i][0].objective = model.reactions[1692]
    oneC_min_models[i][0].reactions[1692].lower_bound = -1000

# To obtain reactions with non-zero fluxes for minimum uptake models
for row in range(len(oneC_min_models)):
    min_model_sols[row].append(oneC_min_models[row][0].optimize())
    for col in range(len(model.reactions)):
        if min_model_sols[row][0].fluxes[col] > 10**-8:
            min_activeRXNs[row][col] = 1

sum_min_activeRXNs = sum(min_activeRXNs)
active_rxnNUMs = []
for i in range(len(sum_min_activeRXNs)):
    if sum_min_activeRXNs[i] == len(oneC_min_models):
        active_rxnNUMs.append(i)


init_core = np.intersect1d(active_rxnNUMs, init_core)
final_core_models = [model.copy() for _ in range(len(init_core) - 1)]
final_core_solutions = []
for i in range(len(init_core) - 1):  # To omit the biomass rxn
    final_core_models[i].reactions[init_core[i]].bounds = (0, 0)
    final_core_solutions.append(final_core_models[i].optimize())
    if solution.fluxes[1692] - 10**-4 <= final_core_solutions[i].fluxes[1692] \
            <= solution.fluxes[1692] + 10**-4:
        init_core = np.delete(init_core, i)
