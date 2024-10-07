import gzip
import os
import numpy as np
from Bio.Seq import Seq
from Bio import Align

# since we have 2 replicates for S108: the first one barcoded with native-barcode01 (file1) and the second with
# native-barcode08(file8). Here, for S108 we will analyze each file independently and combine them together later.
# That's why we have 21 barcoded files (no barcode20) but only 20 unique pMuSICs here.

actual_barcoded_samples = ['barcode01_S108', 'barcode02_S66', 'barcode03_S110', 'barcode04_S92', 'barcode05_S38',
                          'barcode06_S87', 'barcode07_S52', 'barcode08_S108', 'barcode09_S58', 'barcode10_S34',
                          'barcode11_S30', 'barcode12_S56', 'barcode13_S63', 'barcode14_S53', 'barcode15_S47',
                          'barcode16_S71', 'barcode17_S101', 'barcode18_S83', 'barcode19_S59', 'barcode21_S26',
                          'barcode22_S67']

project_root = os.getcwd()
result_path = os.path.join(project_root, 'output/')
os.makedirs(result_path, exist_ok=True)

# step1: decompress .gz files in the fastq_pass package
path = os.path.join(project_root, 'fastq_pass/')

files = os.listdir(path)
files.sort()

original_pool = {}
for file in files:
    if os.path.splitext(file)[1] == ".gz":
        print("\nBelow is the file " + file)
        key = file.split('_')[2]
        print(key)
        with gzip.open(path + file, "rb") as f_in:  # unzip .gz files
            file_content = str(f_in.read())  # change the bytes form into str form
            original_pool[key] = file_content  # add each unzipped .gz file as an element in the pool
            f_in.close()
print("the length of pool_file is " + str(len(original_pool)))
print('the original keys are:')
print(original_pool.keys())

pool = {}
for each_key in original_pool.keys():
    for i in range(len(actual_barcoded_samples)):
        if each_key in actual_barcoded_samples[i]:
            pool[actual_barcoded_samples[i]] = original_pool[each_key]
print('the updated keys of pool are:')
print(pool.keys())

# step2: extract sequence
# Now we have all files unzipped and stored in a dictionary named pool. To analyze each file, consider them as an
# extremely large string, here all sequences are in type "str", not "Seq". Each file contains thousands of reads. Each
# read contains the DNA sequence and other experimental info. Here, we need to extract the sequence of each read.

seq_pool = {}
for key, val in pool.items():
    print(key)
    list_seq_pool = []
    # split the str with fast@420 and return as a list as [junk0, sequence1+junk1, sequence2+junk2...]
    element = val.split("fast@v4.2.0")
    del element[0]  # element[0] is a junk info for the first sequence, now we have the list as [sequence1+junk1,
    # sequence2+junk2...]
    for j in range(len(element)):
        # split the string with the end "+". then for each sub list in pool[i] it looks like [sequence1, junk, ...junk]
        # to extract the sequence, we can only take the first one in list element.
        sequence = element[j].split("+")
        list_seq_pool.append(sequence[0])
    # check if we still have some non-sequence data in the seq_pool
    if "start_time=" in list_seq_pool:
        print("need to del the first element from each file in Pool")
    else:
        print("now you can use seq_pool directly")
    seq_pool[key] = list_seq_pool


# step3: for each fluorescent barcode, filter all sequences within range 2.0kb to 2.8kb into the list_final_pool
final_pool = {}
for key, val in seq_pool.items():
    print(key)
    list_final_pool = []
    for n in range(len(val)):
        count = len(val[n][2:-2])
        if 2000 <= count <= 2800:
            list_final_pool.append(val[n])

    print("The number of extracted sequences is: " + str(len(val)))
    print("The number of filtered sequences is: " + str(len(val)))
    final_pool[key] = list_final_pool

print(final_pool['barcode01_S108'][1])

# now the sequences are ready for analysis.

# step4: Align each read (sequence) to the common cassette (CMV promoter) to re-orientate its direction from 5' to 3',
# start from CMV to TAA

# to reverse the sequence
def rev(it):
    return "".join(reversed(it))


# for each sequence contains 1 position: "CMV-BmaHI-Kozak-fp1"
CMV = "GCTAGCcgttacataacttacggtaaatggcccgcctggctgaccgcccaacgacccccgcccattgacgtcaataatgacgtatgttcccatagtaacgccaataggga\
ctttccattgacgtcaatgggtggagtatttacggtaaactgcccacttggcagtacatcaagtgtatcatatgccaagtacgccccctattgacgtcaatgacggtaaatggcccg\
cctggcattatgcccagtacatgaccttatgggactttcctacttggcagtacatctacgtattagtcatcgctattaccatgGTGATGCGGTTTTGGCAGTACATCAATGGGCGTG\
GATAGCGGTTTGACTCACGGGGATTTCCAAGTCTCCACCCCATTGACGTCAATGGGAGTTTGTTTTGGCACCAAAATCAACGGGACTTTCCAAAATGTCGTAACAACTCCGCCCCAT\
TGACGCAAATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCTACCGCTGATCAGCCTCGTGCttactggcttatcgaaat"

Ori_CMV = Seq(CMV.upper())
C_CMV = Ori_CMV.complement()
Rev_CMV = rev(Ori_CMV)
RC_CMV = Ori_CMV.reverse_complement()

# define a control list for CMV as orientation
orientation_left = [Ori_CMV, C_CMV]
orientation_right = [Rev_CMV, RC_CMV]

aligner = Align.PairwiseAligner()
aligner.mode = "global"

final_seq_pool = {}
final_seq_score = {}
good_alignment_pool_original = {}
for key in sorted(final_pool.keys()):
    print('Below are the scores align with Ori_CMV, C_CMV, Rev_CMV, RC_CMV in final_pool[' + key + ']:')
    print("The number of extracted sequences of " + str(len(seq_pool[key])))
    print("The number of filtered sequences is: " + str(len(final_pool[key])))
    all_scores_cmv = []  # get all 4 scores for each read aligned with CMV
    cmv_score_seq_list = []  # store max score for each seq
    correct_seq_list = []  # collect all re-orientated sequences
    for each_seq in final_pool[key]:
        score_list = []
        seq_left = each_seq[3:800]
        seq_right = each_seq[-800:-3]
        for cmv_l in orientation_left:
            score_left = aligner.score(Seq(seq_left), cmv_l)
            score_list.append(score_left)  # score_list = [score_left] = [Ori_CMV,C_CMV]

        for cmv_r in orientation_right:
            score_right = aligner.score(Seq(seq_right), cmv_r)
            score_list.append(score_right)  # score_list = [score_left + score_right] = [Ori_CMV,C_CMV,Rev_CMV,RC_CMV]

        # The direction of each read could be determined by the max score among 4 values in score_list
        all_scores_cmv.append(score_list)
        max_score = max(score_list)
        cmv_score_seq_list.append(max_score)

        # re-orientate the direction of each read to 5'-3' and store them in the correct_seq_pool
        max_index = np.argmax(score_list)
        each_seq2 = each_seq[3:-3]
        if max_index == 0:
            correct_seq_list.append(each_seq2)
        elif max_index == 1:
            each_sequence = Seq(each_seq2).complement()
            correct_seq_list.append(str(each_sequence))
        elif max_index == 2:
            each_sequence = rev(Seq(each_seq2))
            correct_seq_list.append(str(each_sequence))
        else:
            each_sequence = Seq(each_seq2).reverse_complement()
            correct_seq_list.append(str(each_sequence))

    # a summary of each barcode
    print("All reads aligned with CMV for each MuSIC barcode are " + str(len(final_pool[key])))
    print("the length of cmv_score_seq_list containing the max score of each read is " + str(len(cmv_score_seq_list)))
    # should be the same as the number of the reads of each barcode
    print("The length of correct_seq_list containing filtered sequences (5'-3') with the orientation from CMV to TAA: "
          + str(len(correct_seq_list))) # should be the same as the number of the reads of each barcode
    final_seq_pool.update({key: correct_seq_list})
    final_seq_score.update({key: all_scores_cmv})

    sample_sort = []
    for i in range(len(correct_seq_list)):
        A = [i, cmv_score_seq_list[i], correct_seq_list[i]]
        sample_sort.append(A)
    cmv_score_seq_list_sort = sorted(sample_sort, key=lambda x: x[1], reverse=True)
    print('the length of sorted A is ' + str(len(cmv_score_seq_list_sort))) # should be the same as the number of all
    # reads of this barcode

    # stap5: Classification. Since the length of CMV is 552, any reads over 80% alignment to the CMV region are
    # considered as good reads
    good_alignment = 552 * 0.8
    bad_alignment = 552 * 0.5
    no_alignment = 552 * 0.1

    final_seq_pool_sort = []
    final_seq_pool_1 = []
    final_seq_pool_2 = []
    final_seq_pool_3 = []
    wrong_seq_pool = []
    for each_list in cmv_score_seq_list_sort:
        final_seq_pool_sort.append(each_list[2])
        if each_list[1] >= good_alignment:
            final_seq_pool_1.append(each_list[2])
        elif bad_alignment <= each_list[1] <= good_alignment:
            final_seq_pool_2.append(each_list[2])
        elif no_alignment < each_list[1] <= bad_alignment:
            final_seq_pool_3.append(each_list[2])
        else:
            wrong_seq_pool.append(each_list[2])

    print("the minimal score in cmv_score_seq_pool of this barcode is " + str(min(cmv_score_seq_list)) +
          ', and its index is ' + str(cmv_score_seq_list.index(min(cmv_score_seq_list))))
    print("the maximum score in cmv_score_seq_pool of this barcode is " + str(max(cmv_score_seq_list)) +
          ', and its index is ' + str(cmv_score_seq_list.index(max(cmv_score_seq_list))))

    print("the length of final_seq_pool_1 of " + key + " is " + str(len(final_seq_pool_1)))
    print("the length of final_seq_pool_2 of " + key + " is " + str(len(final_seq_pool_2)))
    print("the length of final_seq_pool_3 of " + key + " is " + str(len(final_seq_pool_3)))
    print("the length of wrong_seq_pool of " + key + " is " + str(len(wrong_seq_pool)))
    good_alignment_pool_original[key] = final_seq_pool_1

# step6 combine the replicates for S108, 'barcode01_S108' and 'barcode08_S108' to prepare the final reads for scoring
good_alignment_pool = {}
merged_list = good_alignment_pool_original['barcode01_S108'] + good_alignment_pool_original['barcode08_S108']
good_alignment_pool['S108'] = merged_list

del good_alignment_pool_original['barcode01_S108']
del good_alignment_pool_original['barcode08_S108']

for key in good_alignment_pool_original.keys():
    new_key = key.split('_')[1]
    print(new_key)
    good_alignment_pool[new_key] = good_alignment_pool_original[key]

print(good_alignment_pool.keys())
print(len(good_alignment_pool))

np.save(result_path + 'Step1_pool.npy', pool, allow_pickle=True)
np.save(result_path + 'Step1_seq_pool.npy', seq_pool, allow_pickle=True)
np.save(result_path + 'Step1_final_pool.npy', final_pool, allow_pickle=True)
np.save(result_path + 'Step1_final_seq_pool.npy', final_seq_pool, allow_pickle=True)
np.save(result_path + 'Step1_good_alignment_pool.npy', good_alignment_pool, allow_pickle=True)


