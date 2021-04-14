import matplotlib.pyplot as plt
import pandas as pd
import os
from pathlib import Path

YEAR = 'YEAR'
CUMULATIVE_CASHFLOW = 'CUMULATIVE_CASHFLOW'

p = str(Path(__file__).parents[3])
cashflow_path = os.path.join(p, r'mfs/processed/econ_run_2021-04-11/base_case_annual_cashflows_payback_period.csv')
df = pd.read_csv(cashflow_path)
fig = plt.figure()

plt.bar(df[YEAR], df[CUMULATIVE_CASHFLOW])
plt.ylabel('Cumulative Cashflow ($M USD)')
plt.title('Cumulative Cashflows Supporting Discounted Payback Period')
plt.xlabel('Year')
save_path = os.path.join(p, r'mfs/processed/econ_run_2021-04-11/payback_period_chart.png')
plt.savefig(save_path)
plt.show()
