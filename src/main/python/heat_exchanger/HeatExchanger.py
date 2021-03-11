import numpy as np


class HeatExchanger:
    """
    For validation of Symmetry heat exchanger
    Shell and Tube for helium and steam heat exchanger
    Use NTU method
    """

    def __init__(self, mass_flow_hot, mass_flow_cold, Cp_hot, Cp_cold, Tin_hot, Tin_cold, U, A, flow_arrangement):
        self.mass_flow_hot = mass_flow_hot
        self.mass_flow_cold = mass_flow_cold
        self.Cp_hot = Cp_hot
        self.Cp_cold = Cp_cold
        self.Tin_hot = Tin_hot
        self.Tin_cold = Tin_cold
        self.U = U
        self.A = A
        self.flow_arr = flow_arrangement

    @staticmethod
    def size_heat_exhanger(heat_duty, heat_transfer_coeff, Tin_hot, Tout_hot, Tin_cold, Tout_cold, tube_OD,
                           tube_pattern='square'):
        pass

    def NTU_method(self):
        C_hot = self.mass_flow_hot * self.Cp_hot
        C_cold = self.mass_flow_cold * self.Cp_cold

        min_C, max_C = min([C_hot, C_cold]), max([C_hot, C_cold])
        rel_C = min_C / max_C

        NTU = (self.U * self.A) / min_C

        # get effectiveness
        epsilon = self.effectiveness_ntu(NTU, rel_C)

        max_heat_flow = min_C * (self.Tin_hot - self.Tin_cold)
        eff_heat_flow = epsilon * max_heat_flow

        Tout_hot = self.Tin_hot - eff_heat_flow / C_hot
        Tout_cold = self.Tin_cold + eff_heat_flow / C_cold

        return {'Tout_hot': Tout_hot,
                'Tout_cold': Tout_cold}

    def effectiveness_ntu(self, NTU, rel_C, N=1):
        """

        Args:
            type:
            rel_C:
            N: <int> number of passes.  Defaults to 1.

        Returns:
            epsilon: <float> effectiveness of heat exchanger
        """

        # https://en.wikipedia.org/wiki/NTU_method
        type_dict = {'Parallel_Flow': lambda NTU, rel_C: (1 - np.exp(-NTU * (1 + rel_C))) / (1 + rel_C),
                     'CounterCurrent_Flow': lambda NTU, rel_C: (1 - np.exp(-NTU * (1 - rel_C))) / (
                             1 - rel_C * np.exp(-NTU * (1 - rel_C)))}

        epsilon = type_dict.get(self.flow_arr, lambda NTU, rel_C: exec('raise(Exception(x))'))(NTU, rel_C)

        # assert that effectiveness is between 0 and 1
        assert 0 <= epsilon <= 1
        return epsilon

    def __call__(self):
        res_dict = self.NTU_method()

        print(res_dict)


# https://www.pdhonline.com/courses/m371/m371content.pdf
def prelim_area_guess(unit_name):
    # todo Joycelyn you need to do something here because we have literature data for our nuclear steam generator

    # https://inis.iaea.org/collection/NCLCollectionStore/_Public/31/057/31057135.pdf?r=1&r=1

    if unit_name == 'SMR_Steam_Generator':
        return 1880  # need to justify assumptions/ limitations around this number
    else:
        raise NotImplemented


if __name__ == '__main__':
    A = prelim_area_guess('SMR_Steam_Generator')

    # INPUTS - UPDATE THESE
    mass_flow_hot = 580  # kg/s
    mass_flow_cold = 130  # kg/s
    Cp_hot = 1.005  # KJ/kg-K
    Cp_cold = 4.18  # KJ/kg-K
    Tin_hot = 523  # K
    Tin_cold = 300  # K
    U = 20  # kW/m2.K
    A = 10  # m2
    flow_arrangement = 'CounterCurrent_Flow'

    x = HeatExchanger(mass_flow_hot, mass_flow_cold, Cp_hot, Cp_cold, Tin_hot, Tin_cold, U, A, flow_arrangement)
    x()

    # todo calculate the heat duty and LMTD for sizing

