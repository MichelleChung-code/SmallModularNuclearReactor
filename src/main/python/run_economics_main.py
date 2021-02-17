from SensitivityAnalysis import SensitivityAnalysis  # in Project settings mark python folder as source
from ProfitabilityAnalysis import ProfitabilityAnalysis
import os
from pathlib import Path
import numpy as np
from datetime import datetime
import econ_constants as const
import pandas as pd


def apply_tolerance(number: object, TOL: object = 1E-10) -> object:
    res = number if abs(number) > TOL else 0
    return res


def process_lower_bound(number):
    number = apply_tolerance(number)

    if number < 0:
        number = 1 + number  # subtracting a negative number

    if number == 0:
        return 1

    return np.floor(number * 10) / 10


def process_upper_bound(number):
    number = apply_tolerance(number)

    if number == 0:
        return 1

    return (np.ceil(number * 10) / 10) + 1


def process_limits(lower_bound, upper_bound):
    limit_ls = [*np.arange(lower_bound, upper_bound + 0.01, 0.1)]

    if not limit_ls:
        return [1]

    return limit_ls


if __name__ == '__main__':
    pd.set_option('display.max_columns', None)

    # inputs
    discount_rate_nominal = 0.15
    inflation_rate = 0.0311
    discount_rate = ((1 + discount_rate_nominal) / (1 + inflation_rate)) - 1

    print('Real Discount Rate: ', discount_rate)

    # read in the base_case revenue and expenses
    p = str(Path(__file__).parents[3])
    base_case_path = os.path.join(p, r'mfs/base_case.csv')
    results_path = os.path.join(p, r'mfs/processed/econ_run_{}'.format(datetime.today().strftime('%Y-%m-%d')))

    if not os.path.exists(results_path):
        os.makedirs(results_path)

    #### BASE CASE INPUTS IN MM USD ####
    TAX_RATE = 0.13
    FCI = 1298.580629
    WC = 210.782767
    LAND = 2.751365
    OFF_SITE_CAPITAL = 64.929031
    START_UP_EXPENSES = 38.957419
    REV = 296.66
    EXPENSES = 65.32

    x = SensitivityAnalysis(base_case_path, results_path, tax_rate=TAX_RATE, FCI=FCI, WC=WC, Land=LAND,
                            i=discount_rate,
                            offsite_capital=OFF_SITE_CAPITAL,
                            start_up_expenses=START_UP_EXPENSES)

    # build base case annual cashflows
    base_case_cashflows = x.build_cashflows()
    base_case_cashflows.to_csv(os.path.join(p, r'mfs/processed/base_case_annual_cashflows.csv'))

    base_case_cashflows = base_case_cashflows['Cash_Flow'].to_frame()
    base_case_cashflows.rename(columns={'Cash_Flow': 'CASHFLOW'}, inplace=True)
    profit_x = ProfitabilityAnalysis(results_path, 'base_case_annual_cashflows_payback_period.csv', base_case_cashflows,
                                     discount_rate)

    net_present_value = profit_x.net_present_value()
    print('NPV: ', net_present_value)

    irr = profit_x.internal_rate_of_return()
    print('IRR: ', irr)

    discounted_payback_period = profit_x.discounted_payback_period()  # defaults to a desired 10% return
    print('Discounted Payback Period', discounted_payback_period)

    # Run Sensitivity Analysis per case
    df = pd.read_csv(os.path.join(p, r'mfs/sensitivity_analysis_cases.csv'), index_col=0)

    dict_map_changes = {const.REVENUE: REV, const.EXPENSES: EXPENSES, const.FCI: FCI}

    for col in const.var_ls:
        df[const.CHANGE_FROM_BASE + '_{}'.format(col)] = round(
            (df[col] - dict_map_changes[col]) / dict_map_changes[col], 2)

    df.to_csv(os.path.join(results_path, 'sensitivity_cases_summary.csv'))

    for case_name in const.CASES_LS:

        index_ls = [case_name + const.LOWER, case_name + const.UPPER]
        new_df = df.loc[index_ls]

        max_change_rev = new_df['CHANGE_FROM_BASE_REVENUE'].max()
        min_change_rev = new_df['CHANGE_FROM_BASE_REVENUE'].min()

        max_change_fci = new_df['CHANGE_FROM_BASE_FCI'].max()
        min_change_fci = new_df['CHANGE_FROM_BASE_FCI'].min()

        max_change_ex = new_df['CHANGE_FROM_BASE_EXPENSES'].max()
        min_change_ex = new_df['CHANGE_FROM_BASE_EXPENSES'].min()

        if any(x < 0 for x in [max_change_rev, max_change_ex, max_change_fci]):
            raise NotImplemented

        # apply tolerance
        max_change_rev = process_upper_bound(max_change_rev)
        min_change_rev = process_lower_bound(min_change_rev)
        max_change_fci = process_upper_bound(max_change_fci)
        min_change_fci = process_lower_bound(min_change_fci)
        max_change_ex = process_upper_bound(max_change_ex)
        min_change_ex = process_lower_bound(min_change_ex)

        adjust_R_LS = process_limits(lower_bound=min_change_rev,
                                     upper_bound=max_change_rev)  # add the 0.01 to the end bound so that it goals to 2.40
        adjust_E_LS = process_limits(lower_bound=min_change_ex,
                                     upper_bound=max_change_ex)  # add the 0.01 to the end bound so that it goals to 1.10
        adjust_FCI_LS = process_limits(lower_bound=min_change_fci,
                                       upper_bound=max_change_fci)  # add the 0.01 to the end bound so that it goals to 1.70

        sensitivity_analysis_results = x(case_name, adjust_R_LS, adjust_E_LS, adjust_FCI_LS)
        sensitivity_analysis_results.to_csv(
            os.path.join(results_path, '{}_sensitivity_analysis_results.csv'.format(case_name)))
