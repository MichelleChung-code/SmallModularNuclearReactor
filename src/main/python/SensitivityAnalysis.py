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

dict_individual_plots = {
    R_AdjustFact: {'title': ' When Varying Annual Revenue', 'x_axis': 'Revenue Adjustment Factor From Base Case',
                   'series_name': 'Varying Annual Revenue', 'legend_loc': 2},
    E_AdjustFact: {'title': ' When Varying Annual Expense', 'x_axis': 'Expense Adjustment Factor From Base Case',
                   'series_name': 'Varying Annual Expense', 'legend_loc': 1},
    FCI_AdjustFact: {'title': ' When Varying Fixed Capital Investment',
                     'x_axis': 'FCI Adjustment Factor From Base Case', 'series_name': 'Varying FCI', 'legend_loc': 1}
}


class SensitivityAnalysis:
    def __init__(self, base_case_path, results_path, tax_rate, FCI, WC, Land, i, offsite_capital, start_up_expenses):
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

    def __call__(self, case_name, adjust_R_LS, adjust_E_LS, adjust_FCI_LS):
        # NOTE:
        # A 0.1 adjustment factor means a 90% decrease from the base case
        # A 1.4 adjustment factor means a 40% increase from the base case
        LS_ALL = [adjust_R_LS, adjust_E_LS, adjust_FCI_LS]

        combinations = list(itertools.product(*LS_ALL))
        combined_results = pd.DataFrame(combinations, columns=['R_AdjustFact', 'E_AdjustFact', 'FCI_AdjustFact'])
        combined_results = self.build_sensitivity_results(combined_results)
        self.plot_combined_results(combined_results, case_name)
        base_case_row = combined_results[
            (abs(combined_results['R_AdjustFact'] - 1) <= 1E-6) & (
                    abs(combined_results['E_AdjustFact'] - 1) <= 1E-6) & (
                    abs(combined_results['FCI_AdjustFact'] - 1) <= 1E-6)]

        self.plot_individual_results(adjust_R_LS=adjust_R_LS, adjust_E_LS=adjust_E_LS, adjust_FCI_LS=adjust_FCI_LS,
                                     base_case_row=base_case_row, case_name=case_name)

        print('==' * 50)
        print('==' * 15, 'CASE: ', case_name.upper(), '==' * 15)
        print('==' * 10, ' COMBINED RESULTS ', '==' * 10)
        assert (combined_results['NPV'].idxmax() == combined_results['IRR'].idxmax())
        assert (combined_results['NPV'].idxmin() == combined_results['IRR'].idxmin())
        print('>>>', ' MAX')
        print(combined_results.loc[combined_results['NPV'].idxmax()])
        print('>>>', ' MIN')
        print(combined_results.loc[combined_results['NPV'].idxmin()])
        print('==' * 50)

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
    def plot_individual_chart(base_case_row, results, output_path, vary_type, case_name):
        """

        :param base_case_row:
        :param results:
        :param output_path:
        :param vary_type: R_AdjustFact, E_AdjustFact, or FCI_AdjustFact
        :return:
        """
        x = results[vary_type]
        y_NPV = results["NPV"]
        y_IRR = results['IRR']

        # Plot NPV
        fig, ax1 = plt.subplots(figsize=(13, 10))

        color = 'tab:red'
        ax1.plot(x, y_NPV, color=color, label=dict_individual_plots[vary_type]['series_name'] + ': NPV')
        base_case_series = ax1.scatter(base_case_row[vary_type], base_case_row['NPV'], c='k', label='Base Case')
        plt.title(dict_individual_plots[vary_type]['title'])
        ax1.set_xlabel(dict_individual_plots[vary_type]['x_axis'])
        ax1.set_ylabel('NPV ($M USD)', color=color)
        ax1.tick_params(axis='y', labelcolor=color)

        ax2 = ax1.twinx()
        color = 'tab:blue'
        ax2.set_ylabel('IRR', color=color)
        ax2.plot(x, y_IRR, color=color, label=dict_individual_plots[vary_type]['series_name'] + ': IRR')
        ax2.scatter(base_case_row[vary_type], base_case_row['IRR'], c='k')
        ax2.tick_params(axis='y', labelcolor=color)

        if vary_type == E_AdjustFact:
            ax1.set_ylim(90, 500)
            ax2.set_ylim(0.120, 0.15)
        elif vary_type == R_AdjustFact:
            ax1.set_ylim(-1000, 3500)
            ax2.set_ylim(0.05, 0.4)
        elif vary_type == FCI_AdjustFact:
            ax1.set_ylim(-700, 800)
            ax2.set_ylim(0.08, 0.24)

        plt.legend([base_case_series], ['Base Case'], loc=dict_individual_plots[vary_type]['legend_loc'])
        plt.tight_layout()
        plt.savefig(os.path.join(output_path, '{}_'.format(case_name) + vary_type + '.png'))

    def helper_populate_remaining_factors_ones(self, results, LS_include):
        LS_fill_ones = [x for x in LS_include if x not in results.columns]
        for col_name in LS_fill_ones:
            results[col_name] = 1

        return results

    def plot_individual_results(self, adjust_R_LS, adjust_E_LS, adjust_FCI_LS, base_case_row, case_name):

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

        results_vary_R.to_csv(
            os.path.join(self.results_path, '{}_sensitivity_analysis_results_vary_R.csv'.format(case_name)))
        results_vary_E.to_csv(
            os.path.join(self.results_path, '{}_sensitivity_analysis_results_vary_E.csv'.format(case_name)))
        results_vary_FCI.to_csv(
            os.path.join(self.results_path, '{}_sensitivity_analysis_results_vary_FCI.csv'.format(case_name)))

        self.plot_individual_chart(base_case_row, results_vary_R, output_path=self.results_path, vary_type=R_AdjustFact,
                                   case_name=case_name)
        self.plot_individual_chart(base_case_row, results_vary_E, output_path=self.results_path, vary_type=E_AdjustFact,
                                   case_name=case_name)
        self.plot_individual_chart(base_case_row, results_vary_FCI, output_path=self.results_path,
                                   vary_type=FCI_AdjustFact, case_name=case_name)

        # print max and min
        LS_DF = [results_vary_R, results_vary_E, results_vary_FCI]
        LS_TYPES = [R_AdjustFact, E_AdjustFact, FCI_AdjustFact]
        for i in range(len(LS_DF)):
            print('==' * 10, LS_TYPES[i], '==' * 10)
            assert (LS_DF[i]['NPV'].idxmax() == LS_DF[i]['IRR'].idxmax())
            assert (LS_DF[i]['NPV'].idxmin() == LS_DF[i]['IRR'].idxmin())
            print('>>>', ' MAX')
            print(LS_DF[i].loc[LS_DF[i]['NPV'].idxmax()])
            print('>>>', ' MIN')
            print(LS_DF[i].loc[LS_DF[i]['NPV'].idxmin()])

    def plot_combined_results(self, results, case_name):
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
        plt.savefig(os.path.join(self.results_path, '{}_NPV_combined.png'.format(case_name)))

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
        plt.savefig(os.path.join(self.results_path, '{}_IRR_combined.png'.format(case_name)))
