import scipy


class Solution:
    """Class for obtaining numeric solutions of models.

    """
    def __init__(self, model, t_0, t_end, step_size):
        """Constructor Method.

        Parameters
        ----------
        model : class
            Class containing equations() function comprising model ODEs
        t_0 : Int
            Start time for simulation
        t_end : Int
            Finish time for simulation
        step_size : Int
            Interval between start and end time at which to return
            simulation results

        """
        self.t_start = t_0
        self.t_end = t_end
        self.t_eval = list(range(t_0, t_end, step_size))
        self.model = model
        if self.model.type == 'no_ab':
            self.equations = model.equations
            self.y0 = [0.36994133339265045, 0,
                       23.41652401591709, 0,
                       1300, 0,
                       0.007671859328948333, 0,
                       0.001518924239131829, 0,
                       0.16065514517469062, 0,
                       0.013527630452127005, 0]
        elif self.model.type == 'no_ab_two_pools':
            self.equations = model.equations
            self.y0 = [1.13, 0.69, 14.4, 0, 1300, 0,
                       5.7e-3, 6.03e-3, 0.952e-3, 0,
                       53.2e-3, 1.21, 2.22e-3, 0]
        elif self.model.type == 'no_ab_silk':
            self.equations = model.equations
            self.y0 = [0.36994133339265045, 0, 23.41652401591709, 0, 1300, 0,
                       0.007671859328948333, 0, 0.001518924239131829, 0,
                       0.16065514517469062, 0, 5.84e-3, 0]
        elif self.model.type == 'one_ab':
            self.equations = model.equations
            self.y0 = [0, 0.36994133339265045, 23.41652401591709, 1300,
                       0, 0.007671859328948333, 0.001518924239131829,
                       0, 0.16065514517469062, 0.013527630452127005,
                       0, 0, 0,
                       0, 0,
                       0, 0]

    def solve(self):
        """Numeric solver.

        Returns
        ----------
        solution: :

        """
        solution = scipy.integrate.solve_ivp(fun=lambda t, y: self.equations(t, y),
                                             t_span=[self.t_eval[0],
                                                     self.t_eval[-1]],
                                             y0=self.y0,
                                             t_eval=self.t_eval,
                                             max_step=1,
                                             method='LSODA')
        return solution
