import pandas as pd
import numpy as np
import numpy_financial as np_fi
from pathlib import Path
import os

# CONSTANTS
CASHFLOW = "CASHFLOW"
DISCOUNTED_CASHFLOW = 'DISCOUNTED_CASHFLOW'
CUMULATIVE_CASHFLOW = 'CUMULATIVE_CASHFLOW'


class ProfitabilityAnalysis:
    def __init__(self, results_path, results_file_name, cashflow_df, discount_rate=0.1):
        self.results_path = results_path
        self.results_file_name = results_file_name
        self.annual_cashflow_df = cashflow_df
        self.annual_discount_rate = discount_rate
        self.annual_cashflow_np_arr = np.array(self.annual_cashflow_df[CASHFLOW])

    def net_present_value(self):
        return round(np_fi.npv(self.annual_discount_rate, self.annual_cashflow_np_arr), 2)

    def internal_rate_of_return(self):
        # discount rate that results in NPV == 0
        return round(np_fi.irr(self.annual_cashflow_np_arr), 4)

    def discounted_payback_period(self):
        # discount future cash flows at desired return rate
        annual_cash_flow_df = self.annual_cashflow_df.copy()
        annual_cash_flow_df[DISCOUNTED_CASHFLOW] = np_fi.pv(self.annual_discount_rate, pmt=0,
                                                            nper=annual_cash_flow_df.index,
                                                            fv=-annual_cash_flow_df[CASHFLOW])
        annual_cash_flow_df[CUMULATIVE_CASHFLOW] = np.cumsum(annual_cash_flow_df[DISCOUNTED_CASHFLOW])

        # get when the project becomes profitable
        mask_positive = annual_cash_flow_df[CUMULATIVE_CASHFLOW] >= 0
        try:
            min_full_year = annual_cash_flow_df[mask_positive].index.values.min() - 1
            profit_point = min_full_year + abs(
                annual_cash_flow_df.loc[min_full_year, CUMULATIVE_CASHFLOW] / annual_cash_flow_df.loc[
                    min_full_year + 1, DISCOUNTED_CASHFLOW])
        except:
            raise Exception("Project not profitable within given time period")

        annual_cash_flow_df.to_csv(os.path.join(self.results_path, self.results_file_name))

        return round(profit_point, 2)

    def EBITDA(self):
        # Need more detailed information
        raise NotImplemented


if __name__ == '__main__':
    print('WARNING - RUNNING FROM CLASS SCRIPT FILE')
    p = str(Path(__file__).parents[3])
    cash_flow_path = os.path.join(p, r'mfs/cashflows_test_mchung.csv')
    cash_flow_df = pd.read_csv(cash_flow_path, index_col=0)
    x = ProfitabilityAnalysis(cash_flow_df)

    net_present_value = x.net_present_value()
    print('NPV: ', net_present_value)

    irr = x.internal_rate_of_return()
    print('IRR: ', irr)

    discounted_payback_period = x.discounted_payback_period()  # defaults to a desired 10% return
    print('Discounted Payback Period', discounted_payback_period)
