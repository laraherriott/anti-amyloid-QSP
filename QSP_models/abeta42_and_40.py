from .parameters import NoAbParameters
# model for two pools of monomers that converge on the same oligomer and plaque pools
# each pool has a different production rate along side different aggregation rates
# for simulating the SILK data these parameters should be the same for both pools, but with one starting from
# steady state (unlabeled) and the other from zero (labeled)
# the production should be the same underlying rate, but the labeled pool is added to at a rate enrichment*synthesis and 
# unlabeled is added to a rate (1-enrichment)*synthesis
# for Abeta 42 and 40 the parameters will be different for the two pools to reflect the different steady state concentrations
# and the different propensities for aggregation


class NoAbModel:
    def __init__(self, leucine=None, Ab40=None, x=None, time_points=None):
        if time_points is None:
            raise ValueError("The total number of timepoints for\
                              the simulation must be provided")
        if x is None:
            self.params = NoAbParameters()
        else:
            self.params = NoAbParameters()
            self.params.k_in = 10**x[0]
            self.params.k_peripheral_production = 10**x[1]
            self.params.k_olig_inc = 10**x[2]
            self.params.k_olig_sep = 10**x[3]
            self.params.k_clear_Abeta_brain = 10**x[4]
            self.params.k_clear_oligomer_brain = 10**x[5]
            self.params.k_monomer_brain_plasma = 10**x[6]
            self.params.k_monomer_plasma_brain = 10**x[7]
            self.params.k_oligomer_brain_plasma = 10**x[8]
            self.params.k_oligomer_plasma_brain = 10**x[9]
            self.params.k_monomer_brain_csf = 10**x[10]
            self.params.k_oligomer_brain_csf = 10**x[11]
            self.params.k_monomer_csf_plasma = 10**x[12]
            self.params.k_oligomer_csf_plasma = 10**x[13]
            self.params.k_clear_Abeta_plasma = 10**x[14]
            self.params.k_clear_oligomer_plasma = 10**x[15]
            self.params.k_monomer_csf_brain = 10**x[16]
            self.params.k_oligomer_csf_brain = 10**x[17]
            self.params.k_olig_inc_ext = 10**x[18]
            self.params.k_olig_sep_ext = 10**x[19]
            self.params.k_plaque_inc = 10**x[20]
            self.params.k_plaque_sep = 10**x[21]
        
        self.type = 'no_ab'
        self.Ab40_agg_scaling = 1
        self.Ab40_sep_scaling = 1
        self.f_leu = [[0]*time_points, [0]*time_points]
        
        if leucine is not None:
            self.type = 'no_ab_silk'
            self.f_leu = leucine  # first element for plasma and so peripheral, second for CSF and so CNS
        
        if Ab40 is not None:
            self.type = 'no_ab_two_pools'
            self.f_leu = [[Ab40[0]]*Ab40[1], [Ab40[0]]*Ab40[1]]  # list of scaling factors for Abeta 40 vs 42 synthesis
            self.Ab40_agg_scaling = Ab40[2]
            self.Ab40_sep_scaling = Ab40[3]

    def equations(self, t, y):
        unlabeled_production = (1-self.f_leu[1][int(t)])*self.params.k_in
        unlabeled_peripheral_production = (1-self.f_leu[0][int(t)])*self.params.k_peripheral_production

        labeled_production = (self.f_leu[1][int(t)])*self.params.k_in
        labeled_peripheral_production = (self.f_leu[0][int(t)])*self.params.k_peripheral_production

        olig_agg_40 = self.params.k_olig_inc * self.Ab40_agg_scaling
        olig_sep_40 = self.params.k_olig_sep * self.Ab40_sep_scaling

        brain_monomer42 = y[0]
        brain_monomer40 = y[1]

        brain_oligomer42 = y[2]
        brain_oligomer40 = y[3]

        brain_plaque42 = y[4]
        brain_plaque40 = y[5]
         
        plasma_monomer42 = y[6]
        plasma_monomer40 = y[7]

        plasma_oligomer42 = y[8]
        plasma_oligomer40 = y[9]

        csf_monomer42 = y[10]
        csf_monomer40 = y[11]
        
        csf_oligomer42 = y[12]
        csf_oligomer40 = y[13]

        d_brain_monomer42 = (unlabeled_production
                             - self.params.k_olig_inc*brain_monomer42
                             + self.params.k_olig_sep*brain_oligomer42
                             - self.params.k_clear_Abeta_brain*brain_monomer42
                             - self.params.k_monomer_brain_plasma*brain_monomer42
                             + self.params.k_monomer_plasma_brain*plasma_monomer42
                             - self.params.k_monomer_brain_csf*brain_monomer42
                             + self.params.k_monomer_csf_brain*csf_monomer42)
        
        d_brain_monomer40 = (labeled_production
                             - olig_agg_40*brain_monomer40
                             + olig_sep_40*brain_oligomer40
                             - self.params.k_clear_Abeta_brain*brain_monomer40
                             - self.params.k_monomer_brain_plasma*brain_monomer40
                             + self.params.k_monomer_plasma_brain*plasma_monomer40
                             - self.params.k_monomer_brain_csf*brain_monomer40
                             + self.params.k_monomer_csf_brain*csf_monomer40)
        
        d_brain_oligomer40 = (olig_agg_40*brain_monomer40
                              - olig_sep_40*brain_oligomer40
                              - self.params.k_plaque_inc*brain_oligomer40
                              + self.params.k_plaque_sep*brain_plaque40
                              - self.params.k_clear_oligomer_brain*brain_oligomer40
                              - self.params.k_oligomer_brain_plasma*brain_oligomer40
                              + self.params.k_oligomer_plasma_brain*plasma_oligomer40
                              - self.params.k_oligomer_brain_csf*brain_oligomer40
                              + self.params.k_oligomer_csf_brain*csf_oligomer40)
        
        d_brain_oligomer42 = (self.params.k_olig_inc*brain_monomer42
                              - self.params.k_olig_sep*brain_oligomer42
                              - self.params.k_plaque_inc*brain_oligomer42
                              + self.params.k_plaque_sep*brain_plaque42
                              - self.params.k_clear_oligomer_brain*brain_oligomer42
                              - self.params.k_oligomer_brain_plasma*brain_oligomer42
                              + self.params.k_oligomer_plasma_brain*plasma_oligomer42
                              - self.params.k_oligomer_brain_csf*brain_oligomer42
                              + self.params.k_oligomer_csf_brain*csf_oligomer42)
        
        d_brain_plaque40 = (self.params.k_plaque_inc*brain_oligomer40
                            - self.params.k_plaque_sep*brain_plaque40)
        d_brain_plaque42 = (self.params.k_plaque_inc*brain_oligomer42
                            - self.params.k_plaque_sep*brain_plaque42)

        d_plasma_monomer42 = (unlabeled_peripheral_production
                              + self.params.k_monomer_brain_plasma*brain_monomer42
                              - self.params.k_monomer_plasma_brain*plasma_monomer42
                              + self.params.k_monomer_csf_plasma*csf_monomer42
                              - self.params.k_monomer_plasma_csf*plasma_monomer42
                              - self.params.k_olig_inc_ext*plasma_monomer42
                              + self.params.k_olig_sep_ext*plasma_oligomer42
                              - self.params.k_clear_Abeta_plasma*plasma_monomer42)

        d_plasma_monomer40 = (labeled_peripheral_production
                              + self.params.k_monomer_brain_plasma*brain_monomer40
                              - self.params.k_monomer_plasma_brain*plasma_monomer40
                              + self.params.k_monomer_csf_plasma*csf_monomer40
                              - self.params.k_monomer_plasma_csf*plasma_monomer40
                              - self.params.k_olig_inc_ext*plasma_monomer40
                              + self.params.k_olig_sep_ext*plasma_oligomer40
                              - self.params.k_clear_Abeta_plasma*plasma_monomer40)

        d_plasma_oligomer40 = (self.params.k_oligomer_brain_plasma*brain_oligomer40
                               - self.params.k_oligomer_plasma_brain*plasma_oligomer40
                               + self.params.k_oligomer_csf_plasma*csf_oligomer40
                               - self.params.k_oligomer_plasma_csf*plasma_oligomer40
                               - self.params.k_clear_oligomer_plasma*plasma_oligomer40)

        d_plasma_oligomer42 = (self.params.k_oligomer_brain_plasma*brain_oligomer42
                               - self.params.k_oligomer_plasma_brain*plasma_oligomer42
                               + self.params.k_oligomer_csf_plasma*csf_oligomer42
                               - self.params.k_oligomer_plasma_csf*plasma_oligomer42
                               - self.params.k_clear_oligomer_plasma*plasma_oligomer42)

        d_csf_monomer42 = (self.params.k_monomer_brain_csf*brain_monomer42
                           - self.params.k_monomer_csf_brain*csf_monomer42
                           - self.params.k_monomer_csf_plasma*csf_monomer42
                           + self.params.k_monomer_plasma_csf*plasma_monomer42
                           - self.params.k_olig_inc_ext*csf_monomer42
                           + self.params.k_olig_sep_ext*csf_oligomer42
                           - self.params.k_clear_Abeta_csf*csf_monomer42)

        d_csf_monomer40 = (self.params.k_monomer_brain_csf*brain_monomer40
                           - self.params.k_monomer_csf_brain*csf_monomer40
                           - self.params.k_monomer_csf_plasma*csf_monomer40
                           + self.params.k_monomer_plasma_csf*plasma_monomer40
                           - self.params.k_olig_inc_ext*csf_monomer40
                           + self.params.k_olig_sep_ext*csf_oligomer40
                           - self.params.k_clear_Abeta_csf*csf_monomer40)

        d_csf_oligomer40 = (self.params.k_oligomer_brain_csf*brain_oligomer40
                            - self.params.k_oligomer_csf_brain*csf_oligomer40
                            - self.params.k_oligomer_csf_plasma*csf_oligomer40
                            + self.params.k_oligomer_plasma_csf*plasma_oligomer40
                            - self.params.k_clear_oligomer_csf*csf_oligomer40)

        d_csf_oligomer42 = (self.params.k_oligomer_brain_csf*brain_oligomer42
                            - self.params.k_oligomer_csf_brain*csf_oligomer42
                            - self.params.k_oligomer_csf_plasma*csf_oligomer42
                            + self.params.k_oligomer_plasma_csf*plasma_oligomer42
                            - self.params.k_clear_oligomer_csf*csf_oligomer42)

        dYdt = [d_brain_monomer42, d_brain_monomer40,
                d_brain_oligomer42, d_brain_oligomer40,
                d_brain_plaque42, d_brain_plaque40,
                d_plasma_monomer42, d_plasma_monomer40,
                d_plasma_oligomer42, d_plasma_oligomer40,
                d_csf_monomer42, d_csf_monomer40,
                d_csf_oligomer42, d_csf_oligomer40]

        return dYdt
