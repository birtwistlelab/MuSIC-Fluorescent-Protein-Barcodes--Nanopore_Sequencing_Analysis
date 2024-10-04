import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
import os

score_of_read = np.load("output/Step2.score_of_read.npy")
print(len(score_of_read))


# step8: create control barcode library for fp combo identification (324), once we have the index of the highest score,
# we know the fp combination

print("Create Control Barcode Library for FP-combo Identification:")
matrix = np.empty((18, 18), dtype=object)
fp_name = ["EBFP2", "mTagBFP2", "mT-Sapphire", "mAmetrine", "mCerulean3", "LSSmOrange", "mBeRFP", "mTFP1", "EGFP",
           "CyOFP1", "mClover3", "mVenus", "mPapaya", "mOrange2", "mRuby3", "mKate2", "mCardinal", "miRFP670"]
counter = 0
ctrl_324 = []
ctrl_barcode = []
for i, fp_m in enumerate(fp_name):
    for j, fp_n in enumerate(fp_name):
        ctrl_324.append([counter, f'{fp_m}', f'{fp_n}'])
        combination = fp_m + '-tPT2A-' + fp_n
        ctrl_barcode.append([counter, combination])
        counter += 1
print(ctrl_barcode)

# step9: identify the indexes of barcodes from the sample pool
score_nested_list = []
for each_list in score_of_read:
    each_score_list = [0] * 324
    max_index = np.argmax(each_list) # identify the MuSIC barcode index of each read
    each_score_list[max_index] = 1
    score_nested_list.append(each_score_list)
score_array = np.array(score_nested_list)
np.save('output/Step3.inferred_score_array_of_pMuSIC_pool.npy', score_array)

# step10: calculate the percentage of each pMuSIC in 324 reference barcode
sum_music = np.sum(score_array, axis=0)  # sum all reads to see the distribution
total_score = np.sum(sum_music)
percentage = sum_music / total_score

print('the maximum fraction of pMuSIC is：')
print(max(percentage))

print("Fraction of each pMuSIC are: ")
percentage.reshape((18, 18))
print(percentage)

data = pd.DataFrame(percentage.reshape((18, 18)))

data.to_excel('output/Step3.pMuSIC_percentage_matrix(18x18).xlsx', index=False)

fp_name = ["EBFP2", "mTagBFP2", "mT-Sapphire", "mAmetrine", "mCerulean3", "LSSmOrange", "mBeRFP", "mTFP1", "EGFP",
           "CyOFP1", "mClover3", "mVenus", "mPapaya", "mOrange2", "mRuby3", "mKate2", "mCardinal", "miRFP670"]

plt.subplots(figsize=(10, 8))
ax = sns.heatmap(data, linewidth=1, square=True, cmap="Purples", xticklabels=fp_name, yticklabels=fp_name,
            linecolor='white', cbar=True, annot=False, vmax=0.015)
plt.subplots_adjust(bottom=0.28)

colorbar = ax.collections[0].colorbar
colorbar.ax.tick_params(labelsize=14)

colorbar.set_label('Fractions of MuSIC Barcodes', fontsize=18, labelpad=12)
ax.tick_params(labelsize=14)
plt.xlabel('Second Fluorescent Protein', fontsize=18, labelpad=12)
plt.ylabel('First Fluorescent Protein', fontsize=18, labelpad=10)


filename = 'output/Step3.fig.2D heatmap for sequencing of pMuSIC pool.png'
plt.savefig(filename, transparent=True)
plt.close()

