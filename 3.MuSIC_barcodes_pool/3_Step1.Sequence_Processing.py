import gzip
import os
import numpy as np
from Bio.Seq import Seq
from Bio import Align

project_root = os.getcwd()
result_path = os.path.join(project_root, 'output/')
os.makedirs(result_path, exist_ok=True)

# step1: decompress .gz files in the fastq_pass package, since it's a mixed MuSIC barcode pool, we will use list instead
# of dictionary here
path = os.path.join(project_root, 'fastq_pass/')
files = os.listdir(path)
files.sort()

pool = []
for file in files:
    if os.path.splitext(file)[1] == ".gz":
        print("\nBelow is the file " + file)
        with gzip.open(path + file, "rb") as f_in:  # unzip .gz files
            file_content = str(f_in.read())  # change the bytes form into str form
            pool.append(file_content)  # add each unzipped .gz file as an element in the pool
            f_in.close()

print("the total number of file is " + str(len(pool)))

# step2: extract sequence
seq_pool = []
for i in range(len(pool)):
    # split the str with fast@420 and return as a list as [junk0, sequence1+junk1, sequence2+junk2...]
    element = pool[i].split("fast@v4.2.0")
    del element[0]  # element[0] is a junk info for the first sequence, now we have the list as [sequence1+junk1,
    # sequence2+junk2...]
    for j in range(len(element)):
        # split the string with the end "+". then for each sub list in pool[i] it looks like [sequence1, junk, ...junk]
        # to extract the sequence, we can only take the first one in list element.
        sequence = element[j].split("+")
        seq_pool.append(sequence[0])

print("The number of extracted sequences is: " + str(len(seq_pool)))
# check if we still have some non-sequence data in the seq_pool
if "start_time=" in seq_pool:
    print("need to del the first element from each file in Pool")
else:
    print("now you can use seq_pool directly")

# step3: for each fluorescent barcode, filter all sequences within range 2.0kb to 2.8kb into the list_final_pool
final_pool = []

for i in range(len(seq_pool)):
    count = len(seq_pool[i])
    if 2000 <= count <= 3000:
        final_pool.append(seq_pool[i])

print("The number of filtered sequences is: " + str(len(final_pool)))
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

all_scores_cmv = []  # get all 4 scores for each read aligned with CMV
cmv_score_seq_pool = [] # store max score for each seq
final_seq_pool = [] # collect all re-orientated sequences

for each_seq in final_pool:  # each sequence waiting for alignment
    score_list = []
    seq_left = each_seq[3:800]
    seq_right = each_seq[-800:-3]
    for cmv_l in orientation_left:
        score_left = aligner.score(Seq(seq_left), cmv_l)
        score_list.append(score_left)

    for cmv_r in orientation_right:
        score_right = aligner.score(Seq(seq_right), cmv_r)
        score_list.append(score_right)

    # The direction of each read could be determined by the max score among 4 values in score_list
    all_scores_cmv.append(score_list)
    max_score = max(score_list)
    cmv_score_seq_pool.append(max_score)

    # re-orientate the direction of each read to 5'-3' and store them in the correct_seq_pool
    max_index = np.argmax(score_list)
    each_seq2 = each_seq[3:-3]
    if max_index == 0:
        final_seq_pool.append(each_seq2)
    elif max_index == 1:
        each_sequence = Seq(each_seq2).complement()
        final_seq_pool.append(str(each_sequence))
    elif max_index == 2:
        each_sequence = rev(Seq(each_seq2))
        final_seq_pool.append(str(each_sequence))
    else:
        each_sequence = Seq(each_seq2).reverse_complement()
        final_seq_pool.append(str(each_sequence))

print("All reads aligned with CMV for each MuSIC barcode are " + str(len(final_pool)))
print("the length of cmv_score_seq_list containing the max score of each read is " + str(len(cmv_score_seq_pool)))
# should be the same as the number of the reads of each barcode
print("The length of correct_seq_list containing filtered sequences (5'-3') with the orientation from CMV to TAA: "
      + str(len(final_seq_pool)))  # should be the same as the number of the reads of each barcode

sample_sort = []
for i in range(len(final_seq_pool)):
    A = [i, cmv_score_seq_pool[i], final_seq_pool[i]]
    sample_sort.append(A)
cmv_score_seq_list_sort = sorted(sample_sort, key=lambda x: x[1], reverse=True)
print('the length of sorted A is ' + str(len(cmv_score_seq_list_sort))) # should be the same as the number of all reads
# of this barcode

# stap5: Classification. Since the length of CMV is 552, any reads over 80% alignment to the CMV region are considered
# as good reads

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


print("the minimal score in cmv_score_seq_pool of this barcode is " + str(min(cmv_score_seq_pool)) +
      ', and its index is ' + str(cmv_score_seq_pool.index(min(cmv_score_seq_pool))))
print("the maximum score in cmv_score_seq_pool of this barcode is " + str(max(cmv_score_seq_pool)) +
      ', and its index is ' + str(cmv_score_seq_pool.index(max(cmv_score_seq_pool))))

print("the length of final_seq_pool_1 is " + str(len(final_seq_pool_1)))
print("the length of final_seq_pool_2 is " + str(len(final_seq_pool_2)))
print("the length of final_seq_pool_3 is " + str(len(final_seq_pool_3)))
print("the length of wrong_seq_pool is " + str(len(wrong_seq_pool)))

np.save(result_path + 'Step1_pool.npy', pool)
np.save(result_path + 'Step1_seq_pool.npy', seq_pool)
np.save(result_path + 'Step1_final_pool.npy', final_pool)
np.save(result_path + 'Step1_final_seq_pool.npy', final_seq_pool)
np.save(result_path + 'Step1_good_alignment_pool.npy', final_seq_pool_1)

