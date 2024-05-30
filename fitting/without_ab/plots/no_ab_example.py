# import necessary packages
import random
import pandas as pd
import matplotlib.pyplot as plt

# import necessary local classes and functions
from QSP_models.abeta42_and_40 import NoAbModel
from QSP_models.solution import Solution

# set seed fpr reproducibility
random.seed(1)

parameters = pd.read_csv('../output/ms_SILK_all_log_no_mon_long.csv', header=None).values
parameters = [x[0] for x in parameters]

model = NoAbModel(x=parameters, time_points=24*364)

# solve the model from 0 to 1 y in time steps of 1 h
solver = Solution(model, 0, int(24*364), 1) 
# will solve with LSODA method
solutions = solver.solve()

time = solutions.t

brain_monomer = solutions.y[0]
brain_oligomer = solutions.y[2]
brain_plaque = solutions.y[4]

plasma_monomer = solutions.y[6]
plasma_oligomer = solutions.y[8]

csf_monomer = solutions.y[10]
csf_oligomer = solutions.y[12]

with open("percentage_changes.txt", "a") as f:
    print('Brain monomer: ', brain_monomer[-1],
          'Brain oligomer: ', brain_oligomer[-1],
          'Brain plaque: ', brain_plaque[-1], file=f)

    print('Brain monomer: ', ((brain_monomer[0]-brain_monomer[-1])/brain_monomer[0])*100, '%', file=f)
    print('Brain oligomer: ', ((brain_oligomer[0]-brain_oligomer[-1])/brain_oligomer[0])*100, '%', file=f)
    print('Brain plaque: ', ((brain_plaque[0]-brain_plaque[-1])/brain_plaque[0])*100, '%', file=f)

    print('Plasma monomer: ', plasma_monomer[-1],
          'Plasma oligomer: ', plasma_oligomer[-1], file=f)

    print('Plasma monomer: ', ((plasma_monomer[0]-plasma_monomer[-1])/plasma_monomer[0])*100, '%', file=f)
    print('Plasma oligomer: ', ((plasma_oligomer[0]-plasma_oligomer[-1])/plasma_oligomer[0])*100, '%', file=f)

    print('CSF monomer: ', csf_monomer[-1],
          'CSF oligomer: ', csf_oligomer[-1], file=f)

    print('CSF monomer: ', ((csf_monomer[0]-csf_monomer[-1])/csf_monomer[0])*100, '%', file=f)
    print('CSF oligomer: ', ((csf_oligomer[0]-csf_oligomer[-1])/csf_oligomer[0])*100, '%', file=f)

fig1 = plt.figure(1)
plt.plot(time, brain_monomer, label='Monomer')
plt.plot(time, brain_oligomer, label='Oligomer')
plt.plot(time, brain_plaque, label='Plaque')
plt.xlabel("Time, months")
plt.ylabel("Species (brain), nM")
plt.legend()
no_months = 12
xticks = [i*(28*24) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
plt.tight_layout()
plt.show()
fig1.savefig('change_in_brain.png', dpi=300)

fig2 = plt.figure(2)
plt.plot(time, plasma_monomer, label='Monomer')
plt.plot(time, plasma_oligomer, label='Oligomer')
plt.xlabel("Time, months")
plt.ylabel("Species (plasma), nM")
plt.legend()
no_months = 12
xticks = [i*(28*24) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
plt.tight_layout()
plt.show()
fig2.savefig('change_in_plasma.png', dpi=300)

fig3 = plt.figure(3)
plt.plot(time, csf_monomer, label='Monomer')
plt.plot(time, csf_oligomer, label='Oligomer')
plt.xlabel("Time, months")
plt.ylabel("Species (CSF), nM")
plt.legend()
no_months = 12
xticks = [i*(28*24) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
plt.tight_layout()
plt.show()
fig3.savefig('change_in_csf.png', dpi=300)
