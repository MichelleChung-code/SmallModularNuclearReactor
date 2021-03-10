class HeatExchanger:
    """
    For validation of Symmetry heat exchanger
    Shell and Tube for helium and steam heat exchanger
    Use NTU method
    """

    def __init__(self, mass_flow_hot, mass_flow_cold, Cp_hot, Cp_cold, Tin_hot, Tin_cold, U, A):
        self.mass_flow_hot = mass_flow_hot
        self.mass_flow_cold = mass_flow_cold
        self.Cp_hot = Cp_hot
        self.Cp_cold = Cp_cold
        self.Tin_hot = Tin_hot
        self.Tin_cold = Tin_cold
        self.U = U
        self.A = A

    @staticmethod
    def size_heat_exhanger(heat_duty, heat_transfer_coeff, Tin_hot, Tout_hot, Tin_cold, Tout_cold, tube_OD,
                           tube_pattern='square'):
        pass

    def NTU_method(self):
        C_hot = self.mass_flow_hot * self.Cp_hot
        C_cold = self.mass_flow_cold * self.Cp_cold

        min_C, max_C = min([C_hot, C_cold]), max([C_hot, C_cold])
        rel_C = min_C / max_C

    def effectiveness_ntu(self):
        pass

    def __call__(self):
        self.NTU_method()


# https://www.pdhonline.com/courses/m371/m371content.pdf
def prelim_area_guess(unit_name):
    # todo Joycelyn you need to do something here because we have literature data for our nuclear steam generator

    # https://inis.iaea.org/collection/NCLCollectionStore/_Public/31/057/31057135.pdf?r=1&r=1

    if unit_name == 'SMR_Steam_Generator':
        return 1880  # need to justify assumptions/ limitations around this number
    else:
        raise NotImplemented


if __name__ == '__main__':
    x = HeatExchanger()
    x()
