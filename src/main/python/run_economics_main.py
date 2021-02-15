from SensitivityAnalysis import SensitivityAnalysis  # in Project settings mark python folder as source
from ProfitabilityAnalysis import ProfitabilityAnalysis
import os
from pathlib import Path
import numpy as np

# inputs
discount_rate_nominal = 0.15
inflation_rate = 0.0311
discount_rate = ((1 + discount_rate_nominal) / (1 + inflation_rate)) - 1

print('Real Discount Rate: ', discount_rate)

# read in the base_case revenue and expenses
p = str(Path(__file__).parents[3])
base_case_path = os.path.join(p, r'mfs/base_case.csv')
results_path = os.path.join(p, r'mfs/processed')

#### BASE CASE INPUTS IN MM USD ####
TAX_RATE = 0.13
FCI = 1298.580629
WC = 210.782767
LAND = 2.751365
OFF_SITE_CAPITAL = 64.929031
START_UP_EXPENSES = 38.957419

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

# Run Sensitivity Analysis
# - H2 price variances

adjust_R_LS = [*np.arange(0.60, 2.41, 0.1)]  # add the 0.01 to the end bound so that it goals to 2.40
adjust_E_LS = [*np.arange(0.90, 1.11, 0.1)]  # add the 0.01 to the end bound so that it goals to 1.10
adjust_FCI_LS = [*np.arange(0.5, 1.71, 0.1)]  # add the 0.01 to the end bound so that it goals to 1.70

sensitivity_analysis_results = x('', adjust_R_LS, adjust_E_LS, adjust_FCI_LS)
sensitivity_analysis_results.to_csv(os.path.join(p, r'mfs/sensitivity_analysis_results.csv'))
