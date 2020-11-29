from SensitivityAnalysis import SensitivityAnalysis  # in Project settings mark python folder as source
from ProfitabilityAnalysis import ProfitabilityAnalysis
import os
from pathlib import Path

# inputs
discount_rate = 0.1

# read in the base_case revenue and expenses
p = str(Path(__file__).parents[3])
base_case_path = os.path.join(p, r'mfs/base_case.csv')
x = SensitivityAnalysis(base_case_path, tax_rate=0.125, FCI=1318.58, WC=214.02, Land=2.75, i=discount_rate,
                        offsite_capital=65.93,
                        start_up_expenses=72.27)

# build base case annual cashflows
base_case_cashflows = x.build_cashflows()
base_case_cashflows.to_csv(os.path.join(p, r'mfs/processed/base_case_annual_cashflows.csv'))

# todo add plotting of the results
base_case_cashflows = base_case_cashflows['Cash_Flow'].to_frame()
base_case_cashflows.rename(columns={'Cash_Flow': 'CASHFLOW'}, inplace=True)
profit_x = ProfitabilityAnalysis(base_case_cashflows, discount_rate)

net_present_value = profit_x.net_present_value()
print('NPV: ', net_present_value)

irr = profit_x.internal_rate_of_return()
print('IRR: ', irr)

discounted_payback_period = profit_x.discounted_payback_period()  # defaults to a desired 10% return
print('Discounted Payback Period', discounted_payback_period)
