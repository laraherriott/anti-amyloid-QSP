# import necessary packages
import random
import pandas as pd
import matplotlib.pyplot as plt

# import necessary local classes and functions
from QSP_models.one_ab_model import OneAbModel
from QSP_models.solution import Solution

# set seed fpr reproducibility
random.seed(1)

parameters = pd.read_csv('../../without_ab/output/ms_SILK_all_log_no_mon_long.csv', header=None).values
parameters = [x[0] for x in parameters]

ab_parameters = pd.read_csv('../output/ms_with_Ab_plasma_csf_log_long_plaque_abs_10.csv', header=None).values
ab_parameters = [x[0] for x in ab_parameters]

model_plaque = OneAbModel(dose=10, bodyweight=83, final_dose=24*364*1.5, dose_interval=14*24,
                          x=parameters, y=ab_parameters)
model_pk = OneAbModel(dose=10, bodyweight=84, final_dose=24*364*1.5, dose_interval=50*24,
                      x=parameters, y=ab_parameters)
model_csf = OneAbModel(dose=15, bodyweight=78, final_dose=24*364*1.5, dose_interval=50*24,
                       x=parameters, y=ab_parameters)

# solve the model from 0 to 1 y in time steps of 1 h
solver_plaque = Solution(model_plaque, 0, int(24*364*1.5), 1)
solver_pk = Solution(model_pk, 0, int(24*30), 1)
solver_csf = Solution(model_csf, 0, int(24*30), 1)

# will solve with LSODA method
solutions_plaque = solver_plaque.solve()
solutions_pk = solver_pk.solve()
solutions_csf = solver_csf.solve()

time_yr = solutions_plaque.t

brain_plaque = [x+y for x, y in zip(solutions_plaque.y[3],
                                    solutions_plaque.y[12])]

plaque_obs = [1300-(0.675*1300), 1300-(0.775*1300)]
time_plaque_obs = [52*7*24, 78*7*24]

with open('lecanemab_plaque_fall.txt', "a") as f:
    print('Change at 52 weeks: ', brain_plaque[52*7*24], '%', file=f)
    print('Change at 78 weeks: ', brain_plaque[-1], '%', file=f)

time_pk_obs = [2*24, 3*24, 11*24, 22*24, 29*24]
pk_obs = [(((((251.9331611)*1e3/1e6)/147181.62))*1e9),
          (((((188.1643137)*1e3/1e6)/147181.62))*1e9),
          (((((49.23056101)*1e3/1e6)/147181.62))*1e9),
          (((((18.76272665)*1e3/1e6)/147181.62))*1e9),
          (((((9.06729014)*1e3/1e6)/147181.62))*1e9)]

time_pk = solutions_pk.t

csf_day1 = (solutions_csf.y[7][24]+solutions_csf.y[15][24]+solutions_csf.y[16][24])
csf_day_10 = (solutions_csf.y[7][24*10]+solutions_csf.y[15][24*10]+solutions_csf.y[16][24*10])
plasma_day_1 = (solutions_csf.y[4][24]+solutions_csf.y[13][24]+solutions_csf.y[14][24])
plasma_day_10 = (solutions_csf.y[4][24*10]+solutions_csf.y[13][24*10]+solutions_csf.y[14][24*10])
csf_prop_day1 = csf_day1/plasma_day_1
csf_prop_day2 = csf_day_10/plasma_day_10

with open('csf_prop.txt', "a") as f:
    print('Ratio at Day 1: ', csf_prop_day1*100, '%', file=f)
    print('Ratio at Day 10: ', csf_prop_day2*100, '%', file=f)

x_error = [24, 24*10]
y_error = [(0.04/100)*plasma_day_1, (0.81/100)*plasma_day_10]

fig1 = plt.figure(1)
plt.plot(time_yr, brain_plaque, label='predicted')
plt.plot(time_plaque_obs, plaque_obs, label='observations', marker='o',
         linestyle="")
plt.xlabel('Time, months')
plt.ylabel('Total plaque concentration, nM')
no_months = 18
xticks = [i*(28*24) for i in range(0, no_months+1, 3)]
xtick_labels = [i for i in range(0, no_months+1, 3)]
plt.xticks(xticks, xtick_labels)
plt.legend()
plt.tight_layout()
plt.show()
fig1.savefig('lecanemab_plaque_fall.png', dpi=300)

fig2 = plt.figure(2)
plt.plot(time_pk, (solutions_pk.y[4]+solutions_pk.y[13]+solutions_pk.y[14]),
         label='predicted')
plt.plot(time_pk_obs, pk_obs, marker='o', label='observed', linestyle="")
plt.xlabel('Time, days')
plt.ylabel('Total antibody concentration in plasma, nM')
no_days = 30
xticks = [i*(24) for i in range(0, no_days+1, 3)]
xtick_labels = [i for i in range(0, no_days+1, 3)]
plt.xticks(xticks, xtick_labels)
plt.show()
fig2.savefig('lecanemab_plasma_pk.png', dpi=300)

fig3 = plt.figure(3)
plt.plot(time_pk, (solutions_csf.y[7]+solutions_csf.y[15]+solutions_csf.y[16]),
         label='predicted')
plt.plot(x_error, y_error, label='observed', marker='o', linestyle="")
plt.xlabel('Time, days')
plt.ylabel('Total antibody concentration in CSF, nM')
no_days = 30
xticks = [i*(24) for i in range(0, no_days+1, 3)]
xtick_labels = [i for i in range(0, no_days+1, 3)]
plt.xticks(xticks, xtick_labels)
plt.show()
fig3.savefig('lecanemab_csf_pk.png', dpi=300)
