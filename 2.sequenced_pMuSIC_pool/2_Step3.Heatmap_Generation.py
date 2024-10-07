import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
import os

project_root = os.getcwd()
result_path = os.path.join(project_root, 'output/')

inferred_barcode_list = np.load(result_path + 'Step2.inferred_barcode_list.npy')
actual_sample_list = np.load(result_path + 'Step2.actual_sample_list.npy')
sample_index_list = np.load(result_path + 'Step2.sample_index_list.npy')
percentage = pd.read_excel(result_path + 'Step2.pMuSIC_percentage_matrix(20x324).xlsx')

# step11: analyze sample fraction
sample_fraction = []
for each_row in percentage.to_numpy():
    barcode_rate_list = []
    for each_idx in sample_index_list:
        barcode_rate = each_row[each_idx]  # get the percentage of each index in the sample_index_list
        barcode_rate_list.append(barcode_rate)
    sample_fraction.append(barcode_rate_list)

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
print(inferred_barcode_list)
print(sample_fraction)
data = pd.DataFrame(sample_fraction)

data.to_excel(result_path + 'Step3.barcode_classification_matrix_for_20x20_heatmap.xlsx')

# now we have the 20x20 sample fraction matrix to classify the sample pMuSICs and the inferred barcodes

rcParams['font.sans-serif'] = ['DejaVu Sans']
rcParams['font.family'] = 'sans-serif'

plt.subplots(figsize=(12, 8))

ax = sns.heatmap(sample_fraction, linewidth=1, square=True, cmap='Purples', xticklabels=inferred_barcode_list,
                 yticklabels=actual_sample_list, linecolor='white', cbar=True, annot=False, vmax=1)
plt.subplots_adjust(left=0.3, bottom=0.3)

colorbar = ax.collections[0].colorbar
colorbar.ax.tick_params(labelsize=14)
ax.tick_params(axis='x', labelsize=10, rotation=90)
plt.xlabel('Inferred', fontsize=18)
plt.ylabel('Actual', fontsize=18)

filename = 'Step3.fig.S6C heatmap for sequencing of 20 sequenced pMuSIC pool.png'
filepath = os.path.join(result_path, filename)
plt.savefig(filepath, transparent=True)
plt.close()


