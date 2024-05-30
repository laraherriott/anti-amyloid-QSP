# import necessary packages
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# import necessary local classes and functions
from QSP_models.abeta42_and_40 import NoAbModel
from QSP_models.solution import Solution

f_leu_plasma = pd.read_csv('../../Data/mean_f_leu.csv', header=None).values
f_leu_plasma = np.append([0], f_leu_plasma)

f_leu_csf = pd.read_csv('../../Data/mean_f_leu_csf.csv', header=None).values
f_leu_csf = np.append([0], f_leu_csf)

#parameters = pd.read_csv('../../MatLab/Fitting/NoAb_fitting/Output_new/SILK_only_all_params_new_MatLab_Apr24.csv', header=None).values
parameters = pd.read_csv('../../MatLab/Fitting/Ab_fitting/ms_SILK_CSF_Plasma100.csv', header=None).values

parameters = [x[0] for x in parameters]

# set seed fpr reproducibility]
random.seed(1)

model = NoAbModel(leucine=[f_leu_plasma, f_leu_csf], x=parameters)

# solve the mode
solver = Solution(model, 0, int(37), 1)
# will solve with LSODA method
solutions = solver.solve()

time = list(solutions.t)

# save unlabeled CSF concetration
csf_monomer = solutions.y[10]
csf_oligomer = solutions.y[12]
# and plasma
plasma_monomer = solutions.y[6]
plasma_oligomer = solutions.y[8]

# save labeled CSF concentration
csf_leu_monomer = solutions.y[11]
csf_leu_oligomer = solutions.y[13]
# and plasma
plasma_leu_monomer = solutions.y[7]
plasma_leu_oligomer = solutions.y[9]

plasma_leu_total = [x+y for x, y in zip(plasma_leu_monomer, plasma_leu_oligomer)]
csf_leu_total = [x+y for x, y in zip(csf_leu_monomer, csf_leu_oligomer)]

# experimentally determined data (mean) to compare to

CSF_times = range(37)

CSF_df = pd.read_csv('../../Data/mean_Ab42_CSF_SILK.csv')
CSF_vals = CSF_df.iloc[:, 0]

error_bars = pd.read_csv('../../Data/silk_ci_std.csv')

x_error = error_bars['time_csf']
std = error_bars['std_csf']
std_adjusted = [((((x*1000)*1e-9)/4514)/1e-9) for x in std]
ci_95 = [2*x for x in std_adjusted]

# plot CSF concentration labeled:
fig1 = plt.figure(1)
plt.plot(time, csf_leu_total, marker='x', label='Predicted')
#plt.plot(CSF_times, CSF_vals, marker='x', label='Observed')
plt.errorbar(x_error, CSF_vals, yerr=ci_95, label='Observed', fmt='o', capsize=5)
plt.xlabel("Time, h")
plt.ylabel("Labeled CSF Amyloid, nM")
plt.legend()
plt.tight_layout()
plt.show()
#fig1.savefig('output/no_ab_new_prefit/CSF_SILK.png', dpi=300)

fig2 = plt.figure(2)
plt.plot(range(37), f_leu_csf[:37], label='CSF')
plt.plot(range(37), f_leu_plasma[:37], label='plasma')
plt.xlabel("Time, h")
plt.ylabel("Labeled leucine fraction")
plt.legend()
plt.show()
#fig2.savefig('output/no_ab_new_prefit/CSF_SILK_leucine.png', dpi=300)
