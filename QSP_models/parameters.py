

class NoAbParameters:
    """Class defining the default parameters for use in class:NoAbModel.

    """
    def __init__(self):
        """Constructor Method.

        """
        self.k_in = 0.026778384278873097
        self.k_peripheral_production = 0.02428272
        self.k_olig_inc = 0.00223718023022605
        self.k_olig_sep = 2.2371802302260502e-06
        self.k_plaque_inc = 0.00010306612852828261
        self.k_plaque_sep = 1.0306612852828261e-10
        self.k_olig_inc_ext = 0.00223718023022605
        self.k_olig_sep_ext = 0.223718023022605
        self.k_clear_Abeta_plasma = 7.100210526315789
        self.k_clear_Abeta_csf = 0
        self.k_clear_Abeta_brain = 0.007163044247787612
        self.k_clear_oligomer_plasma = 0.7085294117647059
        self.k_clear_oligomer_csf = 0
        self.k_clear_oligomer_brain = 2.342083333333333e-05
        self.k_monomer_brain_plasma = 0.004775362831858408
        self.k_monomer_plasma_brain = 1e-10
        self.k_oligomer_brain_plasma = 1.5613888888888887e-05
        self.k_oligomer_plasma_brain = 1e-10
        self.k_monomer_brain_csf = 0.009550725663716815
        self.k_monomer_csf_brain = 1e-10
        self.k_oligomer_brain_csf = 3.1227777777777775e-05
        self.k_oligomer_csf_brain = 1e-10
        self.k_monomer_csf_plasma = 0.077
        self.k_monomer_plasma_csf = 0
        self.k_oligomer_csf_plasma = 0.077
        self.k_oligomer_plasma_csf = 0


class Lecanemab:
    """Class defining the default parameters for use in class:OneAbModel.

    """
    def __init__(self):
        """Constructor Method.

        """
        # antibody specific
        self.onPP = 0.001
        self.offma0 = 2300*self.onPP
        self.offma1 = 0.97*self.onPP
        self.offma2 = 1.8*self.onPP

        # general parameters
        self.k_mAb_plasma_brain = 10**(-9.9999)
        self.k_mAb_brain_plasma = 10**(-4.8067)
        self.k_mAb_brain_csf = 10**(-4.5056)
        self.k_mAb_csf_brain = 10**(-10)
        self.k_mAb_csf_plasma = 10**(-1.1336)
        self.plasma_clearance = 10**(-2.3165)
        self.brain_clearance = 10**(-4.6304)
        self.k_ADCP = 10**(-2.5229)
