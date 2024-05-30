#
# Class to define the model parameters - use as template to generate other parameter files
# note that these parameters are in time units of hour
#

class NoAbParameters:
    def __init__(self):
        self.k_in = 2.71786928e-02  # (initial guess - to be fitted)
        self.k_peripheral_production = 3.05634321e-06 #1.00219530e-02  # (initial guess - to be fitted)
        self.k_olig_inc = 2.65094639e-03  # (initial guess - to be fitted)
        self.k_olig_sep = 1.51906541e-03  # (initial guess - to be fitted)
        self.k_plaque_inc = 0.00010306612852828261
        self.k_plaque_sep = 1.0306612852828261e-10
        self.k_olig_inc_ext = 0
        self.k_olig_sep_ext = 0
        self.k_clear_Abeta_plasma = 7.00151778e+00
        self.k_clear_Abeta_csf = 0
        self.k_clear_Abeta_brain = 6.96300289e-03  # (initial guess - to be fitted)
        self.k_clear_oligomer_plasma = 7.08527642e-01
        self.k_clear_oligomer_csf = 0
        self.k_clear_oligomer_brain = 1.28090355e-04  # (initial guess - to be fitted)
        self.k_monomer_brain_plasma = 2.47169362e-03  # (initial guess - to be fitted)
        self.k_monomer_plasma_brain = 4.88418839e-03  # (initial guess - to be fitted)
        self.k_oligomer_brain_plasma = 2.94675355e-05  # (initial guess - to be fitted)
        self.k_oligomer_plasma_brain = 4.64448574e-08  # (initial guess - to be fitted)
        self.k_monomer_brain_csf = 9.21959269e-03  # (initial guess - to be fitted)
        self.k_monomer_csf_brain = 0
        self.k_oligomer_brain_csf = 2.93657267e-04  # (initial guess - to be fitted)
        self.k_oligomer_csf_brain = 6.09048182e-09
        self.k_monomer_csf_plasma = 7.30564358e-02  # (initial guess - to be fitted)
        self.k_monomer_plasma_csf = 0
        self.k_oligomer_csf_plasma = 7.70076629e-02  # (initial guess - to be fitted)
        self.k_oligomer_plasma_csf = 0


# class NoAbParameters_prefit:
#     def __init__(self):
#         self.k_in = 0.0111  # 0.00456 # (initial guess - to be fitted)
#         self.k_peripheral_production = 0.0061453  # 0.001185 # (initial guess - to be fitted)
#         self.k_olig_inc = 0.001653  # 0.01 # (initial guess - to be fitted)
#         self.k_olig_sep = 0  # 0.00065 # (initial guess - to be fitted)
#         self.k_plaque_inc = 0.000103
#         self.k_plaque_sep = 0
#         self.k_olig_inc_ext = 0
#         self.k_olig_sep_ext = 0
#         self.k_clear_Abeta_plasma = 0.231 # (initial guess - to be fitted)
#         self.k_clear_Abeta_csf = 0
#         self.k_clear_Abeta_brain = 0.00272  # 0.0001168 # (initial guess - to be fitted)
#         self.k_clear_oligomer_plasma = 0.231 # (initial guess - to be fitted)
#         self.k_clear_oligomer_csf = 0
#         self.k_clear_oligomer_brain = 8.9e-6  # 1.53e-5 # (initial guess - to be fitted)
#         self.k_monomer_brain_plasma = 0.00181  # 5.84e-5 # (initial guess - to be fitted)
#         self.k_monomer_plasma_brain = 0  # (initial guess - to be fitted)
#         self.k_oligomer_brain_plasma = 5.94e-6  # 7.64e-6 # (initial guess - to be fitted)
#         self.k_oligomer_plasma_brain = 0  # (initial guess - to be fitted)
#         self.k_monomer_brain_csf = 0.00363  # 5.84e-5 # (initial guess - to be fitted)
#         self.k_monomer_csf_brain = 0
#         self.k_oligomer_brain_csf = 1.19e-5  # 7.64e-6 # (initial guess - to be fitted)
#         self.k_oligomer_csf_brain = 0
#         self.k_monomer_csf_plasma = 0.077  # 0.00124 # (initial guess - to be fitted)
#         self.k_monomer_plasma_csf = 0
#         self.k_oligomer_csf_plasma = 0.077  # 0.0495 # (initial guess - to be fitted)
#         self.k_oligomer_plasma_csf = 0

class NoAbParameters_prefit:
    def __init__(self):
        self.k_in = 0.026778384278873097  # (initial guess - to be fitted)
        self.k_peripheral_production = 0.02428272  # (initial guess - to be fitted)
        self.k_olig_inc = 0.00223718023022605 # (initial guess - to be fitted)
        self.k_olig_sep = 2.2371802302260502e-06 # (initial guess - to be fitted)
        self.k_plaque_inc = 0.00010306612852828261
        self.k_plaque_sep = 1.0306612852828261e-10
        self.k_olig_inc_ext = 0.00223718023022605
        self.k_olig_sep_ext = 0.223718023022605
        self.k_clear_Abeta_plasma = 7.100210526315789 # (initial guess - to be fitted)
        self.k_clear_Abeta_csf = 0
        self.k_clear_Abeta_brain = 0.007163044247787612  # (initial guess - to be fitted)
        self.k_clear_oligomer_plasma = 0.7085294117647059 # (initial guess - to be fitted)
        self.k_clear_oligomer_csf = 0
        self.k_clear_oligomer_brain = 2.342083333333333e-05  #  (initial guess - to be fitted)
        self.k_monomer_brain_plasma = 0.004775362831858408  # (initial guess - to be fitted)
        self.k_monomer_plasma_brain = 1e-10  # (initial guess - to be fitted)
        self.k_oligomer_brain_plasma = 1.5613888888888887e-05  # (initial guess - to be fitted)
        self.k_oligomer_plasma_brain = 1e-10  # (initial guess - to be fitted)
        self.k_monomer_brain_csf = 0.009550725663716815  # (initial guess - to be fitted)
        self.k_monomer_csf_brain = 1e-10
        self.k_oligomer_brain_csf = 3.1227777777777775e-05  # (initial guess - to be fitted)
        self.k_oligomer_csf_brain = 1e-10
        self.k_monomer_csf_plasma = 0.077  # (initial guess - to be fitted)
        self.k_monomer_plasma_csf = 0
        self.k_oligomer_csf_plasma = 0.077  # (initial guess - to be fitted)
        self.k_oligomer_plasma_csf = 0


class OneAbParameters_prefit:
    def __init__(self):
        # antibody specific
        self.onPP = 0.001*360
        self.offma0 = 2290*self.onPP
        self.offma1 = 67.3*self.onPP
        self.offma2 = 1.79*self.onPP

        # general parameters
        self.k_mAb_plasma_brain = 10**(-9.9999)
        self.k_mAb_brain_plasma = 10**(-4.8067)
        self.k_mAb_brain_csf = 10**(-4.5056)
        self.k_mAb_csf_brain = 10**(-10)
        self.k_mAb_csf_plasma = 10**(-1.1336)
        self.plasma_clearance = 10**(-2.3165)
        self.brain_clearance = 10**(-4.6304)
        self.k_ADCP = 10**(-2.5229)


# no antibody parameters fitted with old data and old SILK simulations:
        # self.k_in =  4.50699900e-03 # (initial guess - to be fitted)
        # self.k_peripheral_production = 4.50087350e-04 # (initial guess - to be fitted)
        # self.k_olig_inc = 2.17313773e-02 # (initial guess - to be fitted)
        # self.k_olig_sep = 6.51608755e-04 # (initial guess - to be fitted)
        # self.k_plaque_inc = 0.000103
        # self.k_plaque_sep = 0
        # self.k_olig_inc_ext = 0
        # self.k_olig_sep_ext = 0
        # self.k_clear_Abeta_plasma = 0.231
        # self.k_clear_Abeta_csf = 0
        # self.k_clear_Abeta_brain =  4.74488242e-04 # (initial guess - to be fitted)
        # self.k_clear_oligomer_plasma = 0.231
        # self.k_clear_oligomer_csf = 0
        # self.k_clear_oligomer_brain = 1.68540931e-04 # (initial guess - to be fitted)
        # self.k_monomer_brain_plasma = 2.04675352e-04 # (initial guess - to be fitted)
        # self.k_monomer_plasma_brain = 6.64238956e-07 # (initial guess - to be fitted)
        # self.k_oligomer_brain_plasma = 7.35564688e-06 # (initial guess - to be fitted)
        # self.k_oligomer_plasma_brain = 2.77688030e-03 # (initial guess - to be fitted)
        # self.k_monomer_brain_csf = 2.59947819e-04 # (initial guess - to be fitted)
        # self.k_monomer_csf_brain = 0
        # self.k_oligomer_brain_csf = 7.69621928e-06 # (initial guess - to be fitted)
        # self.k_oligomer_csf_brain = 0
        # self.k_monomer_csf_plasma = 2.09378149e-03 # (initial guess - to be fitted)
        # self.k_monomer_plasma_csf = 0
        # self.k_oligomer_csf_plasma = 4.98133305e-02 # (initial guess - to be fitted)
        # self.k_oligomer_plasma_csf = 0
