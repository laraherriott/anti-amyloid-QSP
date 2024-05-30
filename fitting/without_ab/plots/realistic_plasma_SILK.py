# import necessary packages
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# import necessary local classes and functions
from QSP_models.abeta42_and_40 import NoAbModel
from QSP_models.solution import Solution

f_leu_plasma = [x[0] for x in
                pd.read_csv('../../../../Data/bolus_800_plasma_new.csv',
                            header=None).values]
f_leu_csf = [x[0] for x in
             pd.read_csv('../../../../Data/bolus_800_csf_new.csv',
                         header=None).values]

parameters = pd.read_csv('../output/ms_SILK_all_log_no_mon_long.csv',
                         header=None).values
parameters = [x[0] for x in parameters]

# set seed fpr reproducibility
random.seed(1)

model = NoAbModel(leucine=[f_leu_plasma, f_leu_csf], x=parameters,
                  time_points=25)

# solve the mode
solver = Solution(model, 0, int(25), 1)
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

plasma_leu_total = [x+y for x, y in zip(plasma_leu_monomer,
                                        plasma_leu_oligomer)]

csf_leu_total = [x+y for x, y in zip(csf_leu_monomer,
                                     csf_leu_oligomer)]

plasma_unlabeled_total = [x+y for x, y in zip(plasma_monomer,
                                              plasma_oligomer)]

csf_unlabeled_total = [x+y for x, y in zip(csf_monomer,
                                           csf_oligomer)]

plasma_labeled_fraction = [x/(y+x) for x, y in zip(plasma_leu_total,
                                                   plasma_unlabeled_total)]

csf_labeled_fraction = [x/(y+x) for x, y in zip(csf_leu_total,
                                                csf_unlabeled_total)]

# experimentally determined data (mean) to compare to

plasma_df = pd.read_csv('../../../../Data/plasma_silk_digitized.csv')
plasma_times = plasma_df.iloc[:, 0]
plasma_vals_auc = [(x*0.8799657293261557)*(5.7e-3 + 0.952e-3)
                   for x in plasma_df.iloc[:, 1]]

error_bars = pd.read_csv('../../../../Data/silk_ci_std.csv')
x_error = error_bars['time_plasma']
ci_95 = error_bars[' ci_plasma'][:14]
y_error = plasma_vals_auc[1:]

ci_adjusted = [(x*0.8799657293261557)*(5.7e-3 + 0.952e-3) for x in ci_95]

errors = [y_error[i]-ci_adjusted[i] for i in range(14)]

# plot plasma concentration labeled:
fig1 = plt.figure(1)
plt.plot(time, plasma_leu_total, marker='x', label='Predicted')
plt.errorbar(x_error[:14], y_error, yerr=errors, label='Observed',
             fmt='o', capsize=5)
plt.xlabel("Time, h")
plt.ylabel("Labeled plasma Amyloid, nM")
plt.legend()
plt.tight_layout()
plt.show()
fig1.savefig('Plasma_SILK.png', dpi=300)

fig2 = plt.figure(2)
plt.plot(range(25), f_leu_csf, label='CSF')
plt.plot(range(25), f_leu_plasma, label='plasma')
plt.xlabel("Time, h")
plt.ylabel("Labeled leucine fraction")
plt.legend()
plt.show()
fig2.savefig('Plasma_SILK_leucine.png', dpi=300)

plasma_df_mfl = pd.read_csv('../../../../Data/plasma_silk_digitized.csv')
plasma_vals_auc_mfl = [(x*0.8799657293261557)
                       for x in plasma_df_mfl.iloc[:, 1]]

error_bars_mfl = pd.read_csv('../../../../Data/silk_ci_std.csv')
x_error_mfl = error_bars_mfl['time_plasma']
ci_95_mfl = error_bars_mfl[' ci_plasma'][:14]
y_error_mfl = plasma_vals_auc_mfl[1:]

ci_adjusted_mfl = [(x*0.8799657293261557) for x in ci_95]

errors_mfl = [y_error_mfl[i]-ci_adjusted_mfl[i] for i in range(14)]

fig3 = plt.figure(3)
plt.plot(time, plasma_labeled_fraction, marker='x', label='Predicted')
plt.errorbar(x_error_mfl[:14], y_error_mfl, yerr=errors_mfl, label='Observed',
             fmt='o', capsize=5)
plt.xlabel("Time, h")
plt.ylabel("Plasma Amyloid, MFL")
plt.legend()
plt.tight_layout()
plt.show()
fig3.savefig('Plasma_SILK_MFL.png', dpi=300)
