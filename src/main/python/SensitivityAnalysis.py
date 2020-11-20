import pandas as pd
import numpy as np
from pathlib import Path
import os
import copy

# All values are in $MM (Millions)
R = 'REVENUE'
E = 'EXPENSE'

# Gemma's column headings zzz
FCI = 'FCI'
WC_plus_L = 'WC+L'
d = 'd'  # I think this is depreciation
Profit_b4_tax = 'Profit_before_tax'
Tax = 'Tax'
Cash_flow = 'Cash_Flow'
interest_rate = 'i'


class SensitivityAnalysis:
    def __init__(self, base_case_path, tax_rate=0.5, FCI=10, WC=0.9, Land=0.1, i=0.1):
        self.tax_rate = tax_rate
        self.FCI = FCI
        self.WC = WC
        self.Land = Land
        self.i = i  # interest/discount rate
        self.base_case = pd.read_csv(base_case_path, index_col=0)

    def build_cashflows(self, adjust_R=1, adjust_E=1, adjust_FCI=1):

    def __call__(self):
        pass


if __name__ == '__main__':
    p = str(Path(__file__).parents[3])
    base_case_path = os.path.join(p, r'mfs/base_case.csv')
    x = SensitivityAnalysis(base_case_path)
