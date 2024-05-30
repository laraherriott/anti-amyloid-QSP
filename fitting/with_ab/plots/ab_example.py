# import necessary packages
import random
import pandas as pd
import matplotlib.pyplot as plt

# import necessary local classes and functions
from QSP_models.one_ab_model import OneAbModel
from QSP_models.solution import Solution

# set seed fpr reproducibility]
random.seed(1)

parameters = pd.read_csv('../../MatLab/Fitting/Ab_fitting/ms_SILK_CSF_Plasma100.csv', header=None).values
parameters = [x[0] for x in parameters]

parameters2 = pd.read_csv('../../MatLab/Fitting/Ab_fitting/Output/absplasma_with_Ab_MatLab_May24.csv', header=None).values
parameters2 = [x[0] for x in parameters2]

#model = OneAbModel(dose=10, final_dose=24*364*1.5, dose_interval=14*24, x=parameters, y=parameters2)
model = OneAbModel(dose=10, final_dose=24*364*1.5, dose_interval=50*24, x=parameters, y=parameters2)
#model = OneAbModel(dose=10, final_dose=24*14*6, dose_interval=14*24, x=parameters, y=parameters2)

# solve the model from 0 to 1 y in time steps of 1 h
#solver = Solution(model, 0, int(24*364*1.5), 1) 
solver = Solution(model, 0, int(24*30), 1) 
#solver = Solution(model, 0, int(24*14*8), 1) 

# will solve with LSODA method
solutions = solver.solve()

time = solutions.t

brain_monomer = solutions.y[1]
brain_oligomer = solutions.y[2]
brain_plaque = solutions.y[3]

# print(brain_monomer[-1], brain_oligomer[-1], brain_plaque[-1])

# print('Brain monomer: ', ((brain_monomer[0]-brain_monomer[-1])/brain_monomer[0])*100)
# print('Brain oligomer: ', ((brain_oligomer[0]-brain_oligomer[-1])/brain_oligomer[0])*100)
# print('Brain plaque: ', ((brain_plaque[0]-brain_plaque[-1])/brain_plaque[0])*100)

# plasma_monomer = solutions.y[5]
# plasma_oligomer = solutions.y[6]

# print(plasma_monomer[-1], plasma_oligomer[-1])

# print('Plasma monomer: ', ((plasma_monomer[0]-plasma_monomer[-1])/plasma_monomer[0])*100)
# print('Plasma oligomer: ', ((plasma_oligomer[0]-plasma_oligomer[-1])/plasma_oligomer[0])*100)

# csf_monomer = solutions.y[8]
# csf_oligomer = solutions.y[9]

# print(csf_monomer[-1], csf_oligomer[-1])

# print('CSF monomer: ', ((csf_monomer[0]-csf_monomer[-1])/csf_monomer[0])*100)
# print('CSF oligomer: ', ((csf_oligomer[0]-csf_oligomer[-1])/csf_oligomer[0])*100)

# fig1 = plt.figure(1)
# plt.plot(time, brain_monomer, label='Monomer')
# plt.plot(time, brain_oligomer, label='Oligomer')
# plt.plot(time, brain_plaque, label='Plaque')
# plt.xlabel("Time")
# plt.ylabel("Species (brain), nM")
# plt.legend()
# plt.tight_layout()
# plt.show()
# #fig1.savefig('output/no_ab_new_prefit/SILK_change_in_brain.png', dpi=300)

# fig2 = plt.figure(2)
# plt.plot(time, plasma_monomer, label='Monomer')
# plt.plot(time, plasma_oligomer, label='Oligomer')
# plt.xlabel("Time")
# plt.ylabel("Species (plasma), nM")
# plt.legend()
# plt.tight_layout()
# plt.show()
# #fig2.savefig('output/no_ab_new_prefit/SILK_change_in_plasma.png', dpi=300)

# fig3 = plt.figure(3)
# plt.plot(time, csf_monomer, label='Monomer')
# plt.plot(time, csf_oligomer, label='Oligomer')
# plt.xlabel("Time")
# plt.ylabel("Species (CSF), nM")
# plt.legend()
# plt.tight_layout()
# plt.show()
# #fig3.savefig('output/no_ab_new_prefit/SILK_change_in_csf.png', dpi=300)

# ab_d_obs(1,1) = (((((251.9331611)*1e3/1e6)/147181.62))*1e9);
# ab_d_obs(2,1) = (((((188.1643137)*1e3/1e6)/147181.62))*1e9);
# ab_d_obs(3,1) = (((((49.23056101)*1e3/1e6)/147181.62))*1e9);
# ab_d_obs(4,1) = (((((18.76272665)*1e3/1e6)/147181.62))*1e9);
# ab_d_obs(5,1) = (((((9.06729014)*1e3/1e6)/147181.62))*1e9);


t2 = [2*24, 3*24, 11*24, 22*24, 29*24]
obs = [(((((251.9331611)*1e3/1e6)/147181.62))*1e9), (((((188.1643137)*1e3/1e6)/147181.62))*1e9), (((((49.23056101)*1e3/1e6)/147181.62))*1e9),(((((18.76272665)*1e3/1e6)/147181.62))*1e9),(((((9.06729014)*1e3/1e6)/147181.62))*1e9)  ]


# plasma_ab_error = []
# for i, t in enumerate(t2):
#     error = (solutions.y[4][t] - obs[i])/obs[i]
#     plasma_ab_error.append(error)

# print(plasma_ab_error)
# sq = [x**2 for x in plasma_ab_error]
# print(sum(sq))

# csf_prop_error = []
# for j in range(len(solutions.y[4])-1):
#     error = ((solutions.y[7][j+1]/solutions.y[4][j+1])*100 - 0.1)/0.1
#     csf_prop_error.append(error)

# csf_prop_error = ((solutions.y[7][24+24*14*6]/solutions.y[4][24+24*14*6])*100 - 0.04)
# print((solutions.y[7][24+24*14*6]/solutions.y[4][24+24*14*6])*100, '%')
# print(csf_prop_error)
# #csf_sq = [x**2 for x in csf_prop_error]
# #print(sum(csf_sq**2))
# print(csf_prop_error**2)

# csf_prop_error = ((solutions.y[7][24*10+24*14*6]/solutions.y[4][24*10+24*14*6])*100 - 0.07)
# print((solutions.y[7][24*10+24*14*6]/solutions.y[4][24*10+24*14*6])*100, '%')
# print(csf_prop_error)
# #csf_sq = [x**2 for x in csf_prop_error]
# #print(sum(csf_sq**2))
# print(csf_prop_error**2)

fig4 = plt.figure(4)
plt.plot(time, (solutions.y[4]+solutions.y[13]+solutions.y[14]), label='predicted')
plt.plot(t2, obs, label='observed')
plt.xlabel('Time, h')
plt.ylabel('Antibody concentration, nM')
plt.title('Plasma')
plt.show()

prop1 = (solutions.y[7][24]+solutions.y[15][24]+solutions.y[16][24])/(solutions.y[4][24]+solutions.y[13][24]+solutions.y[14][24])
prop2 = (solutions.y[7][24*10]+solutions.y[15][24*10]+solutions.y[16][24*10])/(solutions.y[4][24*10]+solutions.y[13][24*10]+solutions.y[14][24*10])

fig5 = plt.figure(5)
plt.plot(time, (solutions.y[7]+solutions.y[15]+solutions.y[16]), label='predicted')
plt.xlabel('Time, h')
plt.ylabel('Antibody concentration, nM')
plt.vlines(x=24, ymin=0, ymax=max(solutions.y[7]), color='k', linestyles='dashed')
plt.vlines(x=24*10, ymin=0, ymax=max(solutions.y[7]), color='k', linestyles='dashed')
plt.title('CSF')
plt.show()

fig5 = plt.figure(5)
plt.plot(time, (solutions.y[0]+solutions.y[12]+solutions.y[11]+solutions.y[10]), label='predicted')
plt.xlabel('Time, h')
plt.ylabel('Antibody concentration, nM')
plt.title('Brain')
plt.show()