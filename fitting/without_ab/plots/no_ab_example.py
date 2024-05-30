# import necessary packages
import random
import pandas as pd
import matplotlib.pyplot as plt

# import necessary local classes and functions
from QSP_models.no_ab_model import NoAbModel
from QSP_models.solution import Solution

# set seed fpr reproducibility]
random.seed(1)

#parameters = pd.read_csv('../../MatLab/Fitting/NoAb_fitting/Output_new/SILK_only_all_params_new_MatLab_Apr24.csv', header=None).values
parameters = pd.read_csv('../../MatLab/Fitting/Ab_fitting/ms_SILK_CSF_Plasma100.csv', header=None).values

parameters = [x[0] for x in parameters]

model = NoAbModel(x=parameters)

# solve the model from 0 to 1 y in time steps of 1 h
solver = Solution(model, 0, int(24*364*20), 1) 
# will solve with LSODA method
solutions = solver.solve()

time = solutions.t

brain_monomer = solutions.y[0]
brain_oligomer = solutions.y[1]
brain_plaque = solutions.y[2]

print(brain_monomer[-1], brain_oligomer[-1], brain_plaque[-1])

print('Brain monomer: ', ((brain_monomer[0]-brain_monomer[-1])/brain_monomer[0])*100)
print('Brain oligomer: ', ((brain_oligomer[0]-brain_oligomer[-1])/brain_oligomer[0])*100)
print('Brain plaque: ', ((brain_plaque[0]-brain_plaque[-1])/brain_plaque[0])*100)

plasma_monomer = solutions.y[3]
plasma_oligomer = solutions.y[4]

print(plasma_monomer[-1], plasma_oligomer[-1])

print('Plasma monomer: ', ((plasma_monomer[0]-plasma_monomer[-1])/plasma_monomer[0])*100)
print('Plasma oligomer: ', ((plasma_oligomer[0]-plasma_oligomer[-1])/plasma_oligomer[0])*100)

csf_monomer = solutions.y[5]
csf_oligomer = solutions.y[6]

print(csf_monomer[-1], csf_oligomer[-1])

print('CSF monomer: ', ((csf_monomer[0]-csf_monomer[-1])/csf_monomer[0])*100)
print('CSF oligomer: ', ((csf_oligomer[0]-csf_oligomer[-1])/csf_oligomer[0])*100)

fig1 = plt.figure(1)
plt.plot(time, brain_monomer, label='Monomer')
plt.plot(time, brain_oligomer, label='Oligomer')
plt.plot(time, brain_plaque, label='Plaque')
plt.xlabel("Time")
plt.ylabel("Species (brain), nM")
plt.legend()
plt.tight_layout()
plt.show()
#fig1.savefig('output/no_ab_new_prefit/SILK_change_in_brain.png', dpi=300)

fig2 = plt.figure(2)
plt.plot(time, plasma_monomer, label='Monomer')
plt.plot(time, plasma_oligomer, label='Oligomer')
plt.xlabel("Time")
plt.ylabel("Species (plasma), nM")
plt.legend()
plt.tight_layout()
plt.show()
#fig2.savefig('output/no_ab_new_prefit/SILK_change_in_plasma.png', dpi=300)

fig3 = plt.figure(3)
plt.plot(time, csf_monomer, label='Monomer')
plt.plot(time, csf_oligomer, label='Oligomer')
plt.xlabel("Time")
plt.ylabel("Species (CSF), nM")
plt.legend()
plt.tight_layout()
plt.show()
#fig3.savefig('output/no_ab_new_prefit/SILK_change_in_csf.png', dpi=300)
