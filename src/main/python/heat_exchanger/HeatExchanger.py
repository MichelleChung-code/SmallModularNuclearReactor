import numpy as np
import math


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
    def size_helical_coil_heat_exhanger(N_tubes, shell_Dout, tube_pitch, tube_bundle_height, A, D_tube):
        # "C:\Users\tkdmc\OneDrive - University of Calgary\Capstone_Group25_CHEMENGG\Equipment Spec Sheets\SMR\Literature Sources\Helical_Coil_Steam_Generator_Sizing.pdf"

        shell_Rout = shell_Dout / 2

        shell_Rin = shell_Rout - (N_tubes * tube_pitch / tube_bundle_height)

        # Lengths of the inner, outer, middle layer of tubes
        L_tube_mid = A / (math.pi * N_tubes * D_tube)
        num_rotations = L_tube_mid / (2 * math.pi * (shell_Rin + shell_Rout))
        L_tube_in = 2 * math.pi * shell_Rin * num_rotations
        L_tube_out = 2 * math.pi * shell_Rout * num_rotations

        return {'num_rotations': num_rotations,
                'shell_inner_diameter': shell_Rin * 2}

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

        return res_dict


# https://www.pdhonline.com/courses/m371/m371content.pdf
def prelim_area_guess(unit_name):
    # todo Joycelyn you need to do something here because we have literature data for our nuclear steam generator

    # https://inis.iaea.org/collection/NCLCollectionStore/_Public/31/057/31057135.pdf?r=1&r=1

    if unit_name == 'SMR_Steam_Generator':
        return 1880  # need to justify assumptions/ limitations around this number
    else:
        raise NotImplemented


def log_mean_temperature_difference(flow_type, Tin_hot, Tin_cold, Tout_hot, Tout_cold):
    """
    LMTD = (ΔT1 - ΔT2) / ln(ΔT1/ΔT2)

    Args:
        flow_type: <str> co_current or counter_current

    Returns:
        LMTD: <float> the log mean temperature difference

    """
    if flow_type == 'counter_current':
        delta_T1 = Tin_hot - Tout_cold
        delta_T2 = Tout_hot - Tin_cold
    elif flow_type == 'co_current':
        delta_T1 = Tin_hot - Tin_cold
        delta_T2 = Tout_hot - Tout_cold
    else:
        raise NotImplemented

    return (delta_T1 - delta_T2) / (math.log(delta_T1 / delta_T2))


if __name__ == '__main__':
    A = prelim_area_guess('SMR_Steam_Generator')

    # INPUTS - UPDATE THESE
    mass_flow_hot = 522000 / 3600  # kg/s
    mass_flow_cold = 450880 / 3600  # kg/s
    Cp_hot = (21.064 + 20.918) / 2 / 4  # KJ/kg-K
    Cp_cold = (86.187 + 44.720) / 2 / 18.02  # KJ/kg-K
    Tin_hot = 750 + 273.15  # K
    Tin_cold = 205.3 + 273.15  # K
    U = 1.11e7 / 1000 / A  # kW/m2.K
    flow_arrangement = 'CounterCurrent_Flow'

    x = HeatExchanger(mass_flow_hot, mass_flow_cold, Cp_hot, Cp_cold, Tin_hot, Tin_cold, U, A, flow_arrangement)
    output_temp_dict = x()

    LMTD = log_mean_temperature_difference('counter_current', Tin_hot, Tin_cold, Tout_hot=output_temp_dict['Tout_hot'],
                                           Tout_cold=output_temp_dict['Tout_cold'])
    Q_final = U * A * LMTD

    print('Heat Exchanged (kJ/s): {}'.format(Q_final))

    # sizing
    # todo maybe use this to get the num tubes?
    size_res = HeatExchanger.size_helical_coil_heat_exhanger(N_tubes=182, shell_Dout=2.8, tube_pitch=40e-3,
                                                             tube_bundle_height=10, A=1880, D_tube=25e-3)

    print(size_res)
