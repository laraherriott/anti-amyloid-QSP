#
# Class to solve the differential equations
#
import scipy

class Solution:
    def __init__(self, model, t_0, t_end, step_size):
        self.t_start = t_0
        self.t_end = t_end
        self.t_eval = list(range(t_0, t_end, step_size))
        self.model = model
        if self.model.type == 'no_ab':
            self.equations = model.equations
            #self.y0 = [1.13, 14.4, 1300, 5.7e-3, 0.952e-3, 0.14016, 5.84e-3]
            self.y0 = [0.44645897284219566, 6.951529259897906, 0, 0.0085167614677827, 0.0002597717710968332, 0.16261105599702466, 0.0018235536787239769]
            self.y0 = [0, 0, 0, 0, 0, 0, 0]

        elif self.model.type == 'no_ab_two_pools':
            self.equations = model.equations
            self.y0 = [1.13, 0.69, 14.4, 0, 1300, 0, 5.7e-3, 6.03e-3, 0.952e-3, 0, 53.2e-3, 1.21, 2.22e-3, 0] # for 42/40
        elif self.model.type == 'no_ab_silk':
            self.equations = model.equations
            #self.y0 = [1.13, 0, 14.4, 1300, 5.7e-3, 0, 0.952e-3, 53.2e-3, 0, 2.22e-3] # for SILK
            self.y0 = [1.13, 0, 14.4, 0, 1300, 0, 5.7e-3, 0, 0.952e-3, 0, 0.14016, 0, 5.84e-3, 0]
        elif self.model.type == 'no_ab_silk_sol':
            self.equations = model.equations
            #self.y0 = [1.13, 0, 14.4, 1300, 5.7e-3, 0, 0.952e-3, 53.2e-3, 0, 2.22e-3] # for SILK
            self.y0 = [1.13, 0, 14.4, 1300, (0.14016+5.84e-3), 0, (153.2e-3+2.22e-3), 0]
        elif self.model.type == 'one_ab':
            self.equations = model.equations
            self.y0 = [0, 0.44645897284219566, 6.951529259897906, 1300,
                       0, 0.0085167614677827, 0.0002597717710968332,
                       0, 0.16261105599702466, 0.0018235536787239769,
                       0, 0, 0,
                       0, 0,
                       0, 0]

    def solve(self):
        solution = scipy.integrate.solve_ivp(fun=lambda t, y: self.equations(t, y),
                                             t_span=[self.t_eval[0],
                                                     self.t_eval[-1]],
                                             y0=self.y0,
                                             t_eval=self.t_eval,
                                             max_step=1,
                                             method='LSODA')
        return solution
