import pandas as pd
import numpy as np
from pathlib import Path
import os
import copy
import numpy_financial as np_fi
import itertools
import matplotlib.pyplot as plt

# All values are in $MM (Millions)
R = 'REVENUE'
E = 'EXPENSE'
R_AdjustFact = 'R_AdjustFact'
E_AdjustFact = 'E_AdjustFact'
FCI_AdjustFact = 'FCI_AdjustFact'

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
    def __init__(self, base_case_path, results_path, tax_rate=0.125, FCI=1052.92, WC=170.99, Land=2.75, i=0.1,
                 offsite_capital=52.65,
                 start_up_expenses=31.59):
        self.results_path = results_path
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

        processed_cashflows.loc[1:,
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
        adjust_R_LS = [*np.arange(0.60, 2.4, 0.1)]
        adjust_E_LS = [*np.arange(0.90, 1.1, 0.1)]
        adjust_FCI_LS = [*np.arange(0.5, 1.7, 0.1)]
        LS_ALL = [adjust_R_LS, adjust_E_LS, adjust_FCI_LS]

        combinations = list(itertools.product(*LS_ALL))
        combined_results = pd.DataFrame(combinations, columns=['R_AdjustFact', 'E_AdjustFact', 'FCI_AdjustFact'])
        combined_results = self.build_sensitivity_results(combined_results)
        self.plot_combined_results(combined_results)
        base_case_row = combined_results[
            (combined_results['R_AdjustFact'] == 1) & (combined_results['E_AdjustFact'] == 1) & (combined_results['FCI_AdjustFact'] == 1)]

        self.plot_individual_results(adjust_R_LS=adjust_R_LS, adjust_E_LS=adjust_E_LS, adjust_FCI_LS=adjust_FCI_LS,
                                     base_case_row=base_case_row)

        return combined_results

    def build_sensitivity_results(self, results_df):
        for index, row in results_df.iterrows():
            processed_cashflows = self.build_cashflows(row["R_AdjustFact"], row["E_AdjustFact"], row["FCI_AdjustFact"])
            cashflow_arr = np.array(processed_cashflows[Cash_flow])
            NPV = round(np_fi.npv(self.i, cashflow_arr), 2)
            IRR = round(np_fi.irr(cashflow_arr), 4)
            results_df.loc[index, 'NPV'] = NPV
            results_df.loc[index, 'IRR'] = IRR
        return results_df

    @staticmethod
    def plot_individual_chart(base_case_row, results, type='NPV'):
        pass

    def helper_populate_remaining_factors_ones(self, results, LS_include):
        LS_fill_ones = [x for x in LS_include if x not in results.columns]
        for col_name in LS_fill_ones:
            results[col_name] = 1

        return results

    def plot_individual_results(self, adjust_R_LS, adjust_E_LS, adjust_FCI_LS, base_case_row):

        LS_vary_factors = [R_AdjustFact, E_AdjustFact, FCI_AdjustFact]
        results_vary_R = pd.DataFrame(adjust_R_LS, columns=[R_AdjustFact])
        results_vary_E = pd.DataFrame(adjust_E_LS, columns=[E_AdjustFact])
        results_vary_FCI = pd.DataFrame(adjust_FCI_LS, columns=[FCI_AdjustFact])

        results_vary_R = self.helper_populate_remaining_factors_ones(results_vary_R, LS_vary_factors)
        results_vary_E = self.helper_populate_remaining_factors_ones(results_vary_E, LS_vary_factors)
        results_vary_FCI = self.helper_populate_remaining_factors_ones(results_vary_FCI, LS_vary_factors)

        results_vary_R = self.build_sensitivity_results(results_vary_R)
        results_vary_E = self.build_sensitivity_results(results_vary_E)
        results_vary_FCI = self.build_sensitivity_results(results_vary_FCI)

        results_vary_R.to_csv(os.path.join(self.results_path, 'sensitivity_analysis_results_vary_R.csv'))
        results_vary_E.to_csv(os.path.join(self.results_path, 'sensitivity_analysis_results_vary_E.csv'))
        results_vary_FCI.to_csv(os.path.join(self.results_path, 'sensitivity_analysis_results_vary_FCI.csv'))

    def plot_combined_results(self, results):
        # NPV plot
        fig = plt.figure(figsize=(13, 10))
        ax = fig.add_subplot(111, projection='3d')

        x = np.array(results['R_AdjustFact'])
        y = np.array(results['E_AdjustFact'])
        z = np.array(results['FCI_AdjustFact'])
        c = np.array(results['NPV'])
        ax.set_xlabel('Revenue Fraction of Base Case')
        ax.set_ylabel('Expense Fraction of Base Case')
        ax.set_zlabel('FCI Fraction of Base Case')
        ax.text(1, 1, 1, 'Base Case', color='k')

        plt.title('Net Present Value Sensitivity Analysis Results')

        img = ax.scatter(x, y, z, c=c, cmap='GnBu')
        fig.colorbar(img)
        plt.savefig(os.path.join(self.results_path, 'NPV_combined.png'))
        plt.show()

        fig = plt.figure(figsize=(13, 10))
        ax = fig.add_subplot(111, projection='3d')
        c = np.array(results['IRR'])
        ax.set_xlabel('Revenue Fraction of Base Case')
        ax.set_ylabel('Expense Fraction of Base Case')
        ax.set_zlabel('FCI Fraction of Base Case')
        ax.text(1, 1, 1, 'Base Case', color='k')

        plt.title('Internal Rate of Return Sensitivity Analysis Results')

        img = ax.scatter(x, y, z, c=c, cmap='Reds')
        fig.colorbar(img)
        plt.savefig(os.path.join(self.results_path, 'IRR_combined.png'))
        plt.show()


if __name__ == '__main__':
    print('WARNING - RUNNING FROM CLASS SCRIPT FILE')
    p = str(Path(__file__).parents[3])
    results_path = os.path.join(p, r'mfs/processed')
    base_case_path = os.path.join(p, r'mfs/base_case.csv')
    x = SensitivityAnalysis(base_case_path, results_path)
    results = x()

    results.to_csv(os.path.join(p, r'mfs/sensitivity_analysis_results.csv'))
