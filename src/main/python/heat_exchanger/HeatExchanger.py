import numpy as np
import math
import pprint


class HeatExchanger:
    """
    For validation of Symmetry heat exchanger
    Shell and Tube for helium and steam heat exchanger
    Use NTU method

    Specifically for simplified validation of a Helical Coil Steam Generator
    Symmetry model simplified to a shell and tube heat exchanger
    """

    def __init__(self, mass_flow_hot, mass_flow_cold, Cp_hot, Cp_cold_in, Cp_cold_out, Tin_hot, Tin_cold, U, A,
                 flow_arrangement, h_vap, Tsat, overwrite_heat_bool=False):
        """
        Class definition

        Args:
            mass_flow_hot: <float> mass flow rate of the hot stream
            mass_flow_cold: <float> mass flow rate of the cold stream
            Cp_hot: <float> average mass heat capacity of hot steam
            Cp_cold_in: <float> mass heat capacity of cold stream in
            Cp_cold_out: <float> mass heat capcity of cold stream out
            Tin_hot: <float> temperature of input hot stream, in K
            Tin_cold: <float> temperature of input cold stream, in K
            U: <float> overall heat transfer coefficient W/(m2*K)
            A: <float> heat exchange area m2
            flow_arrangement: <str> flow arrangement type
            h_vap: <float> mass enthalpy of steam vaporization
            Tsat: <float> steam saturation temperature, in K
            overwrite_heat_bool: <float> or <booL> defaults to False, otherwise is a value to use, instead of
            effectiveness-NTU method
        """
        self.mass_flow_hot = mass_flow_hot
        self.mass_flow_cold = mass_flow_cold
        self.Cp_hot = Cp_hot
        self.Cp_cold_in = Cp_cold_in
        self.Cp_cold_out = Cp_cold_out
        self.Tin_hot = Tin_hot
        self.Tin_cold = Tin_cold
        self.U = U
        self.A = A
        self.flow_arr = flow_arrangement
        self.h_vap = h_vap
        self.Tsat = Tsat
        self.overwrite_heat_bool = overwrite_heat_bool

    @staticmethod
    def size_helical_coil_heat_exhanger(rho_shell, mu_shell, Cp_shell, k_shell, mu_tube, Cp_tube, k_tube,
                                        tube_thickness, lambda_type_316, Q, inner_R_shell, LMTD, N_tubes, tube_pitch,
                                        D_tube, shell_mass_flow_rate):
        """
        Run sizing calculations for helical coil steam generator

        Args:
            rho_shell: <float> average mass density of the shell fluid
            mu_shell: <float> average dynamic viscosity of the shell fluid
            Cp_shell: <float> average heat capacity of the shell fluid
            k_shell: <float> average thermal conductivity of the shell fluid
            mu_tube: <float> average dynamic viscosity of the tube fluid
            Cp_tube: <float> average heat capacity of the tube fluid
            k_tube: <float> average thermal conductivity of the tube fluid
            tube_thickness: <float> tube thickness
            lambda_type_316: <float> thermal conductivity of Type-316 stainless steel
            Q: <float> heat transfer flow rate
            inner_R_shell: <float> inner radius of the shell
            LMTD: <float> log mean temperature difference
            N_tubes: <float> number of tubes
            tube_pitch: <float> tube pitch
            D_tube: <float> tube outer diameter
            shell_mass_flow_rate: <float> mass flow rate on the shell side in kg/s

        Returns:
            <dict> containing key sizing results
        """
        # "C:\Users\tkdmc\OneDrive - University of Calgary\Capstone_Group25_CHEMENGG\Equipment Spec Sheets\SMR\Literature Sources\Helical_Coil_Steam_Generator_Sizing.pdf"
        # starting page 28

        # radius of curvature from https://www.researchgate.net/publication/230607485_Investigation_of_Dean_number_and_curvature_ratio_in_a_double-pipe_helical_heat_exchanger
        curve_R = 0.24
        tube_R = D_tube / 2
        rel_t_D = tube_R / curve_R

        # Assume we can scale based on flow rate
        literature_outer_bundle_diameter = 2.8  # m
        flow_rate_scaling_fact = 145 / 96.25
        N_tubes = math.ceil(N_tubes * flow_rate_scaling_fact)

        # todo velocity
        # maximum velocity, shell side
        # use literature diameter for the flow rate cross sectional area
        V_shell_max = shell_mass_flow_rate / (rho_shell * math.pi * (literature_outer_bundle_diameter / 2) ** 2)

        # shell side bulk fluid reynold's number
        Re_b = V_shell_max * D_tube * rho_shell / mu_shell

        # Use largest Nusselt number correlation
        # Assume more than 16 tube rows, such that the correction factor for the number of tube rows approaches 1
        Pr_b = mu_shell * Cp_shell / k_shell
        Nu_b = 0.033 * Re_b ** 0.8 * Pr_b ** 0.36  # neglect the Pr_b/ Pr_w term i.e. don't account for any wall effects

        # calculate convective heat transfer coefficients
        # shell side
        h_shell = k_shell * Nu_b / D_tube  # Using the tube diameter here

        # tube side
        Pr_t = mu_tube * Cp_tube / k_tube
        Nu_straight_t = 0.022 * Re_b ** 0.8 * Pr_t ** 0.5
        # account for helix
        Nu_t = (1 + 3.6 * (1 - rel_t_D) * rel_t_D ** 0.8) * Nu_straight_t
        h_tube = k_tube * Nu_b / D_tube

        # thermal reistance from pipe material: L/lambda
        # L is the thickness of the wall
        # lambda is the thermal conductivity
        # use type-316 stainless steel
        # https://www.tlv.com/global/ME/steam-theory/overall-heat-transfer-coefficient.html
        type_316_resistance = tube_thickness / lambda_type_316

        U = (h_tube ** -1 + type_316_resistance + h_shell ** -1) ** -1

        A = Q / (U * LMTD)

        # Middle layer tube length
        L_tube_mid = A / (math.pi * N_tubes * D_tube)

        # vary the outer shell diameter
        sizing_res_dict = {}

        outer_shell_D_dict = {}
        for outer_shell_diameter_delta in [0.5, 1, 1.5]:
            outer_R_shell = ((inner_R_shell * 2) + outer_shell_diameter_delta) / 2
            str_D = '{}_diameter'.format(outer_R_shell * 2)
            shell_length = tube_pitch * N_tubes / (outer_R_shell - inner_R_shell)
            num_rotations = L_tube_mid / (2 * math.pi * (outer_R_shell + inner_R_shell))
            tube_bundle_height = shell_length / num_rotations
            outer_shell_D_dict[str_D] = {'shell_length': shell_length,
                                         'num_rotations': num_rotations,
                                         'tube_bundle_height': tube_bundle_height}

        sizing_res_dict['shell_outer_diameter'] = outer_shell_D_dict

        sizing_res_dict['A'] = A
        sizing_res_dict['U'] = U
        sizing_res_dict['num_tubes'] = N_tubes

        return sizing_res_dict

    def NTU_method(self):
        """ Runs the effectiveness-NTU calculations """
        C_hot = self.mass_flow_hot * self.Cp_hot
        C_cold = self.mass_flow_cold * self.Cp_cold_in

        min_C, max_C = min([C_hot, C_cold]), max([C_hot, C_cold])
        rel_C = min_C / max_C

        NTU = (self.U * self.A) / min_C

        # get effectiveness
        epsilon = self.effectiveness_ntu(NTU, rel_C)

        max_heat_flow = max_C * (self.Tin_hot - self.Tin_cold)

        # all water converted to steam
        # steam tables https://thermopedia.com/content/1150/
        eff_heat_flow = epsilon * max_heat_flow

        # need to account for latent heat for steam
        # energy balances

        # overwrite if for purely validation purposes
        if self.overwrite_heat_bool:
            eff_heat_flow = self.overwrite_heat_bool

        Tout_hot = self.Tin_hot - eff_heat_flow / C_hot
        Tout_cold = ((eff_heat_flow - (h_vap * self.mass_flow_cold) - (
                self.mass_flow_cold * self.Cp_cold_in * (self.Tsat - self.Tin_cold))) / (
                             self.mass_flow_cold * self.Cp_cold_out)) + self.Tsat

        return {'Tout_hot': Tout_hot,
                'Tout_cold': Tout_cold}

    def effectiveness_ntu(self, NTU, rel_C, N=1):
        """
        Calculates the effectiveness factor in the effectiveness-NTU HE methodology

        Args:
            type: <str> Heat exchanger type
            rel_C: <float> C_min/C_max ratio
            N: <int> number of passes.  Defaults to 1.

        Returns:
            epsilon: <float> effectiveness of heat exchanger
        """

        # https://en.wikipedia.org/wiki/NTU_method
        # https://hyominsite.files.wordpress.com/2015/03/fundamentals-of-heat-and-mass-transfer-6th-edition.pdf section 11.4
        type_dict = {'Parallel_Flow': lambda NTU, rel_C: (1 - np.exp(-NTU * (1 + rel_C))) / (1 + rel_C),
                     'CounterCurrent_Flow': lambda NTU, rel_C: (1 - np.exp(-NTU * (1 - rel_C))) / (
                             1 - rel_C * np.exp(-NTU * (1 - rel_C))),
                     'OneShellPass_Flow': lambda NTU, rel_C: 2 * (1 + rel_C + (1 + (rel_C ** 2) ** 0.5) * (
                             (1 + np.exp(-NTU * (1 + (rel_C ** 2) ** 0.5))) / (
                             1 - np.exp(-NTU * (1 + (rel_C ** 2) ** 0.5))))) ** -1,
                     'CrossSinglePass_Flow': lambda NTU, rel_C: 1 - np.exp(
                         (1 / rel_C) * (NTU ** 0.22) * np.exp(-rel_C * (NTU ** 0.78)) - 1)}

        # CrossSinglePass_Flow is for both fluids unmixed

        epsilon = type_dict.get(self.flow_arr, lambda NTU, rel_C: exec('raise(Exception(x))'))(NTU, rel_C)

        # assert that effectiveness is between 0 and 1
        assert 0 <= epsilon <= 1
        return epsilon

    def __call__(self):
        print('Hot Stream (Helium), Cold Stream (Water/Steam)')
        res_dict = self.NTU_method()

        print(res_dict)
        print({'{} (deg C)'.format(key): value - 273.15 for key, value in res_dict.items()})

        return res_dict


# https://www.pdhonline.com/courses/m371/m371content.pdf
def prelim_area_guess(unit_name):
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
    Q_final = 3.80e5  # kW
    LMTD_sym = log_mean_temperature_difference('counter_current', Tin_hot=750 + 273.15, Tin_cold=205.3 + 273.15,
                                               Tout_hot=250 + 273.15,
                                               Tout_cold=701.5 + 273.15)

    print('UA (W/K): {}'.format(Q_final / LMTD_sym * 1000))
    print('Symmetry LMTD: {}'.format(LMTD_sym))
    A = prelim_area_guess('SMR_Steam_Generator')

    ########################################################################################
    # VALIDATION INPUTS
    mass_flow_hot = 522000 / 3600  # kg/s
    mass_flow_cold = 450880 / 3600  # kg/s
    Cp_hot = (21.064 + 20.920) / 2 / 4  # KJ/kg-K
    Cp_cold_in \
        = 88.173 / 18.02  # KJ/kg-K
    Cp_cold_out = 41.607 / 18.02  # KJ/kg-K
    Tin_hot = 750 + 273.15  # K
    Tin_cold = 205.3 + 273.15  # K
    U = Q_final / LMTD_sym / A  # kW/m2.K
    flow_arrangement = 'CrossSinglePass_Flow'

    # at (2450 + 2550)/2 kPa = 2500 kPa (25 bar)
    h_vap = 1840.2  # kJ/kg from steam tables https://thermopedia.com/content/1150/
    Tsat = 223.989 + 273.15  # K from steam tables

    x = HeatExchanger(mass_flow_hot, mass_flow_cold, Cp_hot, Cp_cold_in, Cp_cold_out, Tin_hot, Tin_cold, U, A,
                      flow_arrangement, h_vap, Tsat, overwrite_heat_bool=Q_final)
    output_temp_dict = x()

    LMTD = log_mean_temperature_difference('counter_current', Tin_hot=Tin_hot, Tin_cold=Tin_cold,
                                           Tout_hot=output_temp_dict['Tout_hot'],
                                           Tout_cold=output_temp_dict['Tout_cold'])
    Q_final = U * A * LMTD

    print('FINAL LMTD, Q:', LMTD, Q_final)

    print('Heat Exchanged (MW): {}'.format(Q_final / 1000))

    ########################################################################################
    # SIZING INPUTS
    # helium is on the shell side
    # water/steam on the tube side
    rho_shell = (3.2545 + 6.1916) / 2  # average density from symmetry in kg/m3
    mu_shell = (4.7035e-5 + 2.8878e-5) / 2  # average dynamic viscosity from symmetry in Pa*s
    Cp_shell = (21.064 + 20.920) / 2 / 4  # heat capacity in kJ/kg*K
    k_shell = (0.3643 + 0.2192) / 2 / 1000  # thermal conductivity kW/mK
    mu_tube = (1.3077e-4 + 3.6669e-5) / 2  # average dynamic viscosity from symmetry in Pa*s
    Cp_tube = (88.238 + 41.511) / 2 / 18.02  # heat capacity in kJ/kg*K
    k_tube = (0.6601 + 0.0948) / 2 / 1000  # thermal conductivity kW/mK
    tube_thickness = 4e-3
    lambda_type_316 = 13 / 1000  # thermal conductivity kW/mK
    inner_R_shell = 2 / 2  # m assume inner diameter of 2
    literature_num_tubes = 182
    tube_outer_diameter = 25e-3
    tube_pitch = 40e-3
    shell_mass_flow_rate = 175  # kg/s

    size_res = HeatExchanger.size_helical_coil_heat_exhanger(rho_shell, mu_shell, Cp_shell, k_shell, mu_tube, Cp_tube,
                                                             k_tube,
                                                             tube_thickness, lambda_type_316, Q_final, inner_R_shell,
                                                             LMTD_sym,
                                                             N_tubes=literature_num_tubes, tube_pitch=tube_pitch,
                                                             D_tube=tube_outer_diameter,
                                                             shell_mass_flow_rate=shell_mass_flow_rate)

    pp = pprint.PrettyPrinter(indent=1)
    pp.pprint(size_res)
