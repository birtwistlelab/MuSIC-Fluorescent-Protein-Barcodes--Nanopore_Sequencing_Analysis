import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
import os

project_root = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join(project_root, 'output/')

inferred_fp_list = np.load(result_path + 'Step2.inferred_fp_list.npy')
actual_fp_list = np.load(result_path + 'Step2.actual_fp_list.npy')
percentage = pd.read_excel(result_path + 'Step2.final_percentage_for_18x18_heatmap.xlsx')

rcParams['font.sans-serif'] = ['DejaVu Sans']
rcParams['font.family'] = 'sans-serif'

plt.subplots(figsize=(8, 6))

ax = sns.heatmap(percentage, linewidth=1, square=True, cmap='Purples', xticklabels=inferred_fp_list,
                 yticklabels=actual_fp_list, linecolor='white', cbar=True, annot=False, vmax=1)
plt.subplots_adjust(bottom=0.28)

colorbar = ax.collections[0].colorbar
colorbar.ax.tick_params(labelsize=14)
plt.xlabel('Inferred', fontsize=18)
plt.ylabel('Actual', fontsize=18)

filename = result_path + 'Step3.fig.S10A heatmap for sequencing of pR-fp pool.png'
plt.savefig(filename, transparent=True)
plt.close()
