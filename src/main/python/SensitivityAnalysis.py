import pandas as pd
import numpy as np
from pathlib import Path
import os
import copy
import numpy_financial as np_fi
import itertools

# All values are in $MM (Millions)
R = 'REVENUE'
E = 'EXPENSE'

# Gemma's column headings zzz
FCI = 'FCI'
WC_plus_L = 'WC+L'
OC_plus_SU = 'OC+SU'  # offsite capital plus start up expenses, not recovered at the end
d = 'd'  # I think this is depreciation
Profit_b4_tax = 'Profit_before_tax'
Tax = 'Tax'
Cash_flow = 'Cash_Flow'
interest_rate = 'i'


class SensitivityAnalysis:
    def __init__(self, base_case_path, tax_rate=0.125, FCI=1052.92, WC=170.99, Land=2.75, i=0.1, offsite_capital=52.65,
                 start_up_expenses=31.59):
        self.tax_rate = tax_rate
        self.FCI = FCI
        self.WC = WC
        self.Land = Land
        self.offsite_cap = offsite_capital
        self.startup_expenses = start_up_expenses
        self.i = i  # interest/discount rate
        self.base_case = pd.read_csv(base_case_path, index_col=0)
        self.base_case_copy = copy.deepcopy(self.base_case)  # self.base_case should not change between scenarios

    def build_cashflows(self, adjust_R=1, adjust_E=1, adjust_FCI=1):
        processed_cashflows = copy.deepcopy(self.base_case)

        processed_cashflows.loc[0, FCI] = -adjust_FCI * self.FCI
        processed_cashflows[R] = self.base_case[R] * adjust_R
        processed_cashflows[E] = self.base_case[E] * adjust_E

        processed_cashflows.loc[0, WC_plus_L] = -(self.Land + self.WC)
        processed_cashflows.loc[max(processed_cashflows.index), WC_plus_L] = (
                self.Land + self.WC)  # assume land and WC recovered at the end

        processed_cashflows.loc[0, OC_plus_SU] = -(self.offsite_cap + self.startup_expenses)

        processed_cashflows[
            d] = adjust_FCI * self.FCI / max(processed_cashflows.index)  # linear relationship for depreciation
        processed_cashflows.loc[1:, Profit_b4_tax] = processed_cashflows[R] + processed_cashflows[E] - \
                                                     processed_cashflows[d]
        processed_cashflows[Tax] = -self.tax_rate * processed_cashflows[Profit_b4_tax]

        processed_cashflows[Cash_flow] = processed_cashflows[[FCI, R, E, WC_plus_L, OC_plus_SU, Tax]].sum(axis=1)

        return processed_cashflows.fillna(0)

    def __call__(self):
        # NOTE:
        # A 0.1 adjustment factor means a 90% decrease from the base case
        # A 1.4 adjustment factor means a 40% increase from the base case
        adjust_R_LS = [*np.arange(0.1, 2.1, 0.1)]  # ranging from a 90% decrease to a doubling
        adjust_E_LS = copy.deepcopy(adjust_R_LS)
        adjust_FCI_LS = copy.deepcopy(adjust_R_LS)
        LS_ALL = [adjust_R_LS, adjust_E_LS, adjust_FCI_LS]

        combinations = list(itertools.product(*LS_ALL))
        results = pd.DataFrame(combinations, columns=['R_AdjustFact', 'E_AdjustFact', 'FCI_AdjustFact'])
        for index, row in results.iterrows():
            processed_cashflows = self.build_cashflows(row["R_AdjustFact"], row["E_AdjustFact"], row["FCI_AdjustFact"])
            cashflow_arr = np.array(processed_cashflows[Cash_flow])
            NPV = round(np_fi.npv(self.i, cashflow_arr), 2)
            IRR = round(np_fi.irr(cashflow_arr), 4)
            results.loc[index, 'NPV'] = NPV
            results.loc[index, 'IRR'] = IRR

        return results


if __name__ == '__main__':
    p = str(Path(__file__).parents[3])
    base_case_path = os.path.join(p, r'mfs/base_case.csv')
    x = SensitivityAnalysis(base_case_path)
    results = x()

    results.to_csv(os.path.join(p, r'mfs/sensitivity_analysis_results.csv'))

# todo add plotting of the results
