import pandas as pd
import os
from pathlib import Path


def initial_values_node_num(initial_vals_max_node_df, node_num):
    """
    Args:
        initial_vals_max_node_df: <pd.DataFrame> of the Initial Values for 10 Nodes
        node_num: <int> number of nodes to keep the initial value for

    Returns:
        df <pd.DataFrame> filtered dataframe only containing the values up to the desired node number
    """
    if node_num > 10:
        raise Exception('Max node number is 10')

    # Max will be 10
    if node_num == 10:
        return initial_vals_max_node_df

    node_str_ls = ['Node {}'.format(node_num+1)]

    for i in range(node_num + 2, 11):
        node_str_ls.append('Node {}'.format(i))

    df = initial_vals_max_node_df[~initial_vals_max_node_df[1].isin(node_str_ls)]

    return df

if __name__ == '__main__':
    p = os.path.join(str(Path(__file__).parents[1]), 'matlab')
    p_file = os.path.join(p, 'Initial Values_10_Nodes.csv')

    initial_vals_max_node_df = pd.read_csv(p_file, header=None)

    for N in range(1, 10):
        filtered_df = initial_values_node_num(initial_vals_max_node_df, node_num=N)
        p_result = os.path.join(p, 'Initial Values_{}_Nodes.csv'.format(N))
        filtered_df.to_csv(p_result, index=False, header=False)


