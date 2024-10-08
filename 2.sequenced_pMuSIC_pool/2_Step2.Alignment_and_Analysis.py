import numpy as np
from Bio import Align
from Bio.Seq import Seq
import pandas as pd
import time
import multiprocessing as mp
import os

project_root = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join(project_root, 'output/')

# step6: create the control sequence library for alignment
fp_1 = "ATGGTgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgaggggcgagggcgagggcgatgccaccaacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgagccacggcgtgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcacctacaagacccgcgccgaggtgaagttcgagggcgacaccctagtgaaccgcatcgagctgaagggcgtcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaacttcaacagccacaacatctatatcatggccgtcaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacgtggaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacagccactacctgagcacccagtccgtgctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttccgcaccgccgccgggatcactctcgGCATGGACGAGCTGTACAAG"
fp_2 = "ATGGTgtctaagggcgaagagctgattaaggagaacatgcacatgaagctgtacatggagggcaccgtggacaaccatcacttcaagtgcacatccgagggcgaaggcaagccctacgagggcacccagaccatgagaatcaaggtggtcgagggcggccctctccccttcgccttcgacatcctggctactagcttcctctacggcagcaagaccttcatcaaccacacccagggcatccccgacttcttcaagcagtccttccctgagggcttcacatgggagagagtcaccacatacgaagatgggggcgtgctgaccgctacccaggacaccagcctccaggacggctgcctcatctacaacgtcaagatcagaggggtgaacttcacatccaacggccctgtgatgcagaagaaaacactcggctgggaggccttcaccgagacactgtaccccgctgacggcggcctggaaggcagaaacgacatggccctgaagctcgtgggcgggagccatctgatcgcaaacgccaagaccacatatagatccaagaaacccgctaagaacctcaagatgcctggcgtctactatgtggactacagactggaaagaatcaaggaggccaacaacgagacatacgtcgagcagcacgaggtggcagtggccagatactgcgacctccctagcaaactggggcacaagcttaat"
fp_3 = "ATGGTgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccttctcctacggcgtgatggtgttcgcccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaacttcaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggccaacttcaagatccgccacaacatcgaggacggcggcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcatccagtccaagctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctggGCATGGACGAGCTGTACAAG"
fp_4 = "ATGGTgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgcgcggcgagggcgagggcgatgccaccaacggcaagctgaccctgaagttcatctgcacctccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgtcttacggcgtgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatctccttcaaggacgacggcagctacaggacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaacatgaacgtgtgggacgcgtatatcacggccgacaagcagaagaacggcatcaaagcgaacttcaagatcgagcacaacgtcgaggacggcggcgtgcagctcgccgacgcgtaccagcagaacacccccatcggcgacggctccgtgctgctgcctgacaaccactacctgagcttccagagcaagctgttcaaagaccccaacgagcagcgcgatcacatggtcctgctggagttcgttaccgccgccgggatcactcccgGCATGGACGAGCTGTACAAG"
fp_5 = "ATGGTGAGCAAGGGCGAGGAGctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgagctggggcgtgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaacgccatccacggcaacgtctatatcaccgccgacaagcagaagaacggcatcaaggccaacttcggcctcaactgcaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccaagctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcgGCATGGACGAGCTGTACAAG"
fp_6 = "ATGGTgagcaagggcgaggagaataacatggccatcatcaaggagttcatgcgcttcaaggtgcgcatggagggctccgtgaacggccacgagttcgagatcgagggcgagggcgagggccgcccctacgagggctttcagaccgttaagctgaaggttaccaagggtggaccactgcctttcgcctgggacatattgtcacctcagttcacctacggctccaaggcctacgtgaagcaccccgccgacatccccgactacctcaagctgtccttccccgagggcttcaagtgggagcgcgtgatgaacttcgaggatggcggcgtggtgaccgtgactcaggactcctccctgcaggacggcgagttcatctacaaggtgaagctgcgcggcaccaacttcccctccgacggacccgtaatgcagaagaaaaccatgggcatggaagcctcatctgaacggatgtaccccgaggacggcgcactgaagggcgaggacaagctcaggctgaagctgaaggacggcggccactacacctccgaggtcaagaccacctacaaggccaagaagcccgtgcagttgccaggcgcctacatcgtcgacatcaagttggacatcacctcacacaacgaggactacaccatcgtggaacagtacgaacgcgccgagggccgccactccaccggcgGCATGGACGAGCTGTACAAG"
fp_7 = "ATGGTgtctaagggcgaagagctgattaaggagaacatgcatatgaagctgtacatggagggcaccgtgaacaaccaccacttcaagtgcacatccgagggcgaaggcaagccctacgagggcacccagaccatgagaatcaaggtggtcgagggcggccctctccccttcgccttcgacatcctggctaccagcttcatgtacggcagcaaaaccttcatcaaccacacccagggcatccccgacttctttaagcagtccttccctgagggcttcacatgggagagatccaccacatacgaggacgggggcgtgctgaccgctacccaggacaccagcctccaggacggctgcctcatctacaacgtcaagatcagaggggtgaacttcccatccaacggccctgtgatgcagaagaaaacactcggttgggaagcctccacagagatgctgtacccagctgacggaggtctggaaggaagagactacatggctctgaagctcgtgggcgggggccacctgatctgcaacgccaagaccacatacagatccaagaaacccgctaagaacctcaagatgcccggcgtctactatgtggacaggagactggaaagaatcaaggaggccgacaaagagacaagcgtcgagcagcacgaggtggctgtggccagatactgcgacctccctagcaaactggggcacaaa"
fp_8 = "atggtgagcaagggcgaggagactacaatgggcgtaatcaagcccgacatgaagatcaagctgaagatggagggcaacgtgaatggccacgccttcgtgatcgagggcgagggcgagggcaagccctacgacggcaccaacaccatcaacctggaggtgaaggagggagcccccctgcccttctcctacgacattctgaccaccgcgttcgcctacggcaacagggccttcaccaagtaccccgacgacatccccaactacttcaagcagtccttccccgagggctactcttgggagcgcaccatgaccttcgaggacaagggcatcgtgaaggtgaagtccgacatctccatggaggaggactccttcatctacgagatacacctcaagggcgagaacttcccccccaacggccccgtgatgcagaagaaaaccaccggctgggacgcctccaccgagaggatgtacgtgcgcgacggcgtgctgaagggcgacgtcaagcacaagctgctgctggagggcggcggccaccaccgcgttgacttcaagaccatctacagggccaagaaggcggtgaagctgcccgactatcactttgtggaccaccgcatcgagatcctgaaccacgacaaggactacaacaaggtgaccgtttacgagagcgccgtggcccgcaactccaccgacggcatggacgagctgtacaag"
fp_9 = "atggtgagcaagggcgaggaGctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccaagctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcgGCATGGACGAGCTGTACAAG"
fp_10 = "ATGGTgagcaagggcgaggagctgatcaaggagaacatgagaagcaagctgtacctggaaggcagcgtgaacggccaccagttcaagtgcacccacgaaggggagggcaagccctacgagggcaagcagaccaacaggatcaaggtggtggagggaggccccctgccgttcgcattcgacatcctggccacccactttatgtacgggagcaaggtgttcatcaagtaccccgccgacctccccgattattttaagcagtccttccctgagggcttcacatgggagagagtcatggtgttcgaggacgggggcgtgctgaccgccacccaggacaccagcctccaggacggcgagctcatctacaacgtcaaggtcagaggggtgaacttcccagccaacggccccgtgatgcagaagaaaacactgggctgggaaccaagcaccgaaaccatgtacccagctgacggaggtctggaaggaagatgcgacaaggccctgaagctcgtgggaggtggacacttacacgtcaacttcaagaccacatacaagtccaagaaacccgtgaagatgcccggcgtccactacgtggaccgcagactggaaagaatcaaggaggccgacaacgagacttacgtcgagcagtacgagcacgctgtggccagatactccaacctgggtggagGCATGGACGAGCTGTACAAG"
fp_11 = "atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtccgcggcgagggcgagggcgatgccaccaacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccttcggctacggcgtggcctgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatctctttcaaggacgacggtacctacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaacttcaacagccactacgtctatatcacggccgacaagcagaagaactgcatcaaggctaacttcaagatccgccacaacgttgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagccatcagtccaagctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggattacacatgGCATGGACGAGCTGTACAAG"
fp_12 = "atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagctgatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgggctacggcctgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcaccgccgacaagcagaagaacggcatcaaggccaacttcaagatccgccacaacatcgaggacggcggcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagctaccagtccaagctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctctacaag"
fp_13 = "ATGGTgagcaagggcgaggggCAATCGAAGCATGGACTAAAGGAGGAGATGACGGTGAAGTACCACATGGAAGGTTGCGTTAATGGTCACAAATTTGTCATTACTGGAGAGGGCATTGGAAACCCTTTTAAGGGTAAGCAAACCGCAAATTTGTGTGTGATAGAAGGAGGCCCGCTGCCGTTCTCGGAGGACATTTTAAGTCCGGGTTTTAAATATGGTGACCGGATTTTCACAGAGTACCCACAGGATATTGTAGATTACTTCAAGAACTCATGTCCAGCGGGCTATACGTGGGAAAGGAGCTACCTCTTTGAGGACGGAGCGGTCTGTCGATGCAACGTGGACATAACAGTCTCTGAAAAGGAGAACTGCATCTATCACAAAAGTATCTTCAGAGGGGTGAATTTTCCCGCCGACGGCCCCGTAATGAAGAAAATGACCACTAATTGGGAAGCTAGTACCGAGAAAATTGTGCCTGTTCCAAAGCAAGGGATATTAAAGGGAAAGGTCAAAATGTGCCTGTTGCTGAAGGATGGCGGTCGTTATCATTGCCAGTTTGATACGGTATATAAAGCTAAATCAGTGCCCTCCAAGATGCCAGAATGGCATTTTATACAGCATAAGCTCCTTAGGGAGGACCGCTCTGATGCTAAAAACCAGAAATGGCAATTAACAGAACATGCAATCGCAGGCATGGACGAGCTGTACAAG"
fp_14 = "ATGGTgagcaagggcgaggagaataacatggccatcatcaaggagttcatgcgcttcaaggtgcgcatggagggctccgtgaacggccacgagttcgagatcgagggcgagggcgagggccgcccctacgagggctttcagaccgctaagctgaaggtgaccaagggtggccccctgcccttcgcctgggacatcctgtcccctcatttcacctacggctccaaggcctacgtgaagcaccccgccgacatccccgactacttcaagctgtccttccccgagggcttcaagtgggagcgcgtgatgaactacgaggacggcggcgtggtgaccgtgacccaggactcctccctgcaggacggcgagttcatctacaaggtgaagctgcgcggcaccaacttcccctccgacggccccgtgatgcagaagaaaaccatgggctgggaggcctcctccgagcggatgtaccccgaggacggtgccctgaagggcaagatcaagatgaggctgaagctgaaggacggcggccactacacctccgaggtcaagaccacctacaaggccaagaagcccgtgcagctgcccggcgcctacatcgtcgacatcaagttggacatcacctcccacaacgaggactacaccatcgtggaacagtacgaacgcgccgagggccgccactccaccggcgGCATGGACGAGCTGTACAAG"
fp_15 = "atggtgtctaagggcgaagagctgatcaaggaaaatatgcgtatgaaggtggtcatggaaggttcggtcaacggccaccaattcaaatgcacaggtgaaggagaaggcagaccgtacgagggagtgcaaaccatgaggatcaaagtcatcgagggaggacccctgccatttgcctttgacattcttgccacgtcgttcatgtatggcagccgtacctttatcaagtacccggccgacatccctgatttctttaaacagtcctttcctgagggttttacttgggaaagagttacgagatacgaagatggtggagtcgtcaccgtcacgcaggacaccagccttgaggatggcgagctcgtctacaacgtcaaggtcagaggggtaaactttccctccaatggtcccgtgatgcagaagaaaaccaagggttgggagcctaatacagagatgatgtatccagcagatggtggtctgagaggatacactgacatcgcactgaaagttgatggtggtggccatctgcactgcaacttcgtgacaacttacaggtcaaaaaagaccgtcgggaacatcaagatgcccggtgtccatgccgttgatcaccgcctggaaaggatcgaggagagtgacaatgaaacctacgtagtgcaaagagaagtggcagttgccaaatacagcaaccttggtggtgGCATGGACGAGCTGTACAAG"
fp_16 = "ATGGTgagcgagctgattaaggagaacatgcacatgaagctgtacatggagggcaccgtgaacaaccaccacttcaagtgcacatccgagggcgaaggcaagccctacgagggcacccagaccatgagaatcaaggcggtcgagggcggccctctccccttcgccttcgacatcctggctaccagcttcatgtacggcagcaaaaccttcatcaaccacacccagggcatccccgacttctttaagcagtccttccccgagggcttcacatgggagagagtcaccacatacgaagatgggggcgtgctgaccgctacccaggacaccagcctccaggacggctgcctcatctacaacgtcaagatcagaggggtgaacttcccatccaacggccctgtgatgcagaagaaaacactcggctgggaggcctccaccgagacactgtaccccgctgacggcggcctggaaggcagagccgacatggccctgaagctcgtgggcgggggccacctgatctgcaacttaaagaccacatacagatccaagaaacccgctaagaacctcaagatgcccggcgtctactatgtggacagaaggctggaaagaatcaaggaggccgacaaagagacttacgtcgagcagcacgaggtggctgtggccagatactgcgacctccctagcaaactggggcacaga"
fp_17 = "ATGGTgagcaagggcgaggagctgatcaaggagaacatgcacatgaagctgtacatggaaggcaccgtgaacaaccaccacttcaagtgcaccaccgaaggggagggcaagccctacgagggcacccagacccagaggattaaggtggtggagggaggccccctgccgttcgcattcgacatcctggccacctgctttatgtacgggagcaagaccttcatcaaccacacccagggcatccccgatttctttaagcagtccttccctgagggcttcacatgggagagagtcaccacatacgaggacgggggcgtgcttaccgttacccaggacaccagcctccaggacggctgcttgatctacaacgtcaagctcagaggggtgaacttcccatccaacggccctgtgatgcagaagaaaacactcggctgggaggccaccaccgagactctgtaccccgctgacggcggcctggaaggcagatgcgacatggccctgaagctcgtgggcgggggccacctgcactgcaaccttaagaccacatacagatccaagaaacccgctaagaacctcaagatgcccggcgtctactttgtggaccgcagactggaaagaatcaaggaggccgacaatgagacatacgtcgagcagcacgaggtggctgtggccagatactgcgacctccctagcaaactggggcacaaacttaatgGCATGGACGAGCTGTACAAG"
fp_18 = "ATGGTagcaggtcatgcctctggcagccccgcattcgggaccgcctctcattcgaattgcgaacatgaagagatccacctcgccggctcgatccagccgcatggcgcgcttctggtcgtcagcgaacatgatcatcgcgtcatccaggccagcgccaacgccgcggaatttctgaatctcggaagcgtactcggcgttccgctcgccgagatcgacggcgatctgttgatcaagatcctgccgcatctcgatcccaccgccgaaggcatgccggtcgcggtgcgctgccggatcggcaatccctctacggagtactgcggtctgatgcatcggcctccggaaggcgggctgatcatcgaactcgaacgtgccggcccgtcgatcgatctgtcaggcacgctggcgccggcgctggagcggatccgcacggcgggttcactgcgcgcgctgtgcgatgacaccgtgctgctgtttcagcagtgcaccggctacgaccgggtgatggtgtatcgtttcgatgagcaaggccacggcctggtattctccgagtgccatgtgcctgggctcgaatcctatttcggcaaccgctatccgtcgtcgactgtcccgcagatggcgcggcagctgtacgtgcggcagcgcgtccgcgtgctggtcgacgtcacctatcagccggtgccgctggagccgcggctgtcgccgctgaccgggcgcgatctcgacatgtcgggctgcttcctgcgctcgatgtcgccgtgccatctgcagttcctgaaggacatgggcgtgcgcgccaccctggcggtgtcgctggtggtcggcggcaagctgtggggcctggttgtctgtcaccattatctgccgcgcttcatccgtttcgagctgcgggcgatctgcaaacggctcgccgaaaggatcgcgacgcggatcaccgcgcttgagagc"

NheI_cmv_kozak = "CTAGCcgttacataacttacggtaaatggcccgcctggctgaccgcccaacgacccccgcccattgacgtcaataatgacgtatgttcccatagtaacgccaatagggactttccattgacgtcaatgggtggagtatttacggtaaactgcccacttggcagtacatcaagtgtatcatatgccaagtacgccccctattgacgtcaatgacggtaaatggcccgcctggcattatgcccagtacatgaccttatgggactttcctacttggcagtacatctacgtattagtcatcgctattaccatgGTGATGCGGTTTTGGCAGTACATCAATGGGCGTGGATAGCGGTTTGACTCACGGGGATTTCCAAGTCTCCACCCCATTGACGTCAATGGGAGTTTGTTTTGGCACCAAAATCAACGGGACTTTCCAAAATGTCGTAACAACTCCGCCCCATTGACGCAAATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCTACCGCTGATCAGCCTCGTGCttactggcttatcgaaatgGATCCGCCGCCACC"
tPT2A_kozak = "GGcTCCggcgccacaaacttctctctgctaaagcaagcaggtgatgttgaagaaaaccccgggcctggcagcggcgagggcaggggaagtctacttacatgcggggacgtggaggaaaatcccggcccaGGATCCGCCGCCACC"
tail_MfpI = "TAAACTGTATGATTCTAGATAATGaattctgaaacataaaatgaatgcaatt" # optional

fp = []
for i in range(1, 19):
    name = "fp_"+str(i)
    fp.append(eval(name))  # change the str "fp_1" to variable fp_1
    fp[i-1] = fp[i-1].upper()

cmv_fp = []  # list "cmv_fp" contains NheI-CMV-BamHI-Kozak-ATG-fp
for i in range(len(fp)):
    sequence = NheI_cmv_kozak.upper() + fp[i]
    cmv_fp.append(sequence)

cmv_music = [] # list "cmv_music" contains NheI-cmv-BamHI-kozak-fpm-tPT2A-fpn-MfeI
for each_seq in cmv_fp:
    for i in range(len(fp)):
        music_seq = each_seq + tPT2A_kozak.upper() + fp[i]
        cmv_music.append(music_seq)

cmv_music = np.array(cmv_music)
np.save(result_path + 'Step2.cmv_music.npy', cmv_music)

# step7: align and score to cmv_music to minimize the global scoring errors
final_good_read_pool = np.load(result_path + 'Step1_good_alignment_pool.npy', allow_pickle=True).item()

aligner = Align.PairwiseAligner()
aligner.mode = "global"


def scoring(task_dict):
    result = {}
    for key, value in task_dict.items():
        print("start to analyze pMuSIC sample: " + key)
        start_time = time.time()
        readable_start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))
        print(f"start_time of {key}: {readable_start_time}")
        music_score_list = []
        print("the length of total sequences in pMuSIC_" + key + " is ", len(value))
        for each_seq in value:
            score_list = []
            for each_ctrl in cmv_music:
                each_score = aligner.score(Seq(each_seq[:-2]), Seq(each_ctrl))
                score_list.append(each_score)
            music_score_list.append(score_list)
        music_score_array = np.array(music_score_list)
        result.update({key: music_score_array})
        print('the scoring results for ' + key + " are:")
        print('the size of music_score_array is ', music_score_array.size)
        print('the shape of music_score_array is ', music_score_array.shape)
        end_time = time.time()
        readable_end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))
        print(f"end_time of {key}: {readable_end_time}")
    return result


if __name__ == '__main__':
    print("Test the Control Sequence Library for Alignment:")
    print("the length of fp1 is " + str(len(fp[0])))
    print("the sequence of fp1 is:")
    print(fp[0])
    print("the length of list fp is " + str(len(fp)))
    print("the length of NheI-cmv-BamHI-kozak-fpm is:")
    print(len(cmv_fp[0]))
    print("the length of NheI-cmv-BamHI-kozak-fpm-tPT2A-fpn-MfeI is:")
    print(len(cmv_music[0]))
    print("the length of the control sequence library is:")
    print(len(cmv_music))

    print("Align Each pMuSIC Read in Each Sample Group to the Control Sequence Library:")

    start = time.perf_counter()
    num_cores = int(mp.cpu_count()) # num_cores = 8

    #  divide 20 sequenced pMuSICs into num_cores (4) of tasks, each contains 20/num_cores (5) pMuSICs.
    task_list = [{k: final_good_read_pool[k] for k in list(final_good_read_pool.keys())[i:i+5]}
                 for i in range(0, len(final_good_read_pool), 5)]

    with mp.Pool(processes=num_cores) as pool:
        results = pool.map(scoring, task_list)

    # Merge the results into the final result dictionary: score_of_read
    score_of_read = {}
    for result in results:
        score_of_read.update(result)

    print("final result of scoring:")
    for key, value in score_of_read.items():
        print(f"{key}: array of shape {value.shape}")

    np.save(result_path + "Step2.score_of_read.npy", score_of_read, allow_pickle=True)
    print("now you have score_of_read stored!")

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
    actual_sample_list = [
        "mBeRFP-tPT2A-mCerulean3(S108)", "EBFP2-tPT2A-mOrange2(S66)", "mT-Sapphire-tPT2A-mCerulean3(S110)",
        "mTagBFP2-tPT2A-mKate2(S92)", "mAmetrine-tPT2A-LSSmOrange(S38)", "mCerulean3-tPT2A-mTFP1(S87)",
        "LSSmOrange-tPT2A-mBeRFP(S52)", "mBeRFP-tPT2A-mTFP1(S58)", "mTFP1-tPT2A-mCerulean3(S34)",
        "EGFP-tPT2A-mTagBFP2(S30)", "EGFP-tPT2A-mCerulean3(S56)", "mClover3-tPT2A-mTFP1(S63)",
        "CyOFP1-tPT2A-mCardinal(S53)",
        "mVenus-tPT2A-mCerulean3(S47)", "mPapaya-tPT2A-EBFP2(S71)", "mOrange2-tPT2A-mAmetrine(S101)",
        "mRuby3-tPT2A-mCerulean3(S83)", "mRuby3-tPT2A-CyOFP1(S59)", "mCardinal-tPT2A-mKate2(S26)",
        "miRFP670-tPT2A-mClover3(S67)"
    ]
    np.save(result_path + 'Step2.actual_sample_list.npy', actual_sample_list)

    sample_index_list = []
    barcode_info_list = []
    for i in range(len(actual_sample_list)):
        for idx, each_barcode in ctrl_barcode:
            if each_barcode in actual_sample_list[i]:
                sample_index_list.append(idx)
                barcode_info_list.append([idx, each_barcode])
    print("indentify the index list of actual sample barcodes in control barcode library:")
    print(sample_index_list)
    print("to verify the index list of actual sample barcodes:")
    print('barcode_info_list', barcode_info_list)
    print('actual_sample_list', actual_sample_list)
    np.save(result_path + 'Step2.sample_index_list.npy', sample_index_list)
    np.save(result_path + 'Step2.sample_barcode_info_list.npy', barcode_info_list)

    score_of_reads_nested_list = []
    inferred_barcode_list = []
    for key in list(score_of_read.keys()):
        print(key)
        print(len(score_of_read[key]))
        score_array = score_of_read[key]
        score_list = score_array.tolist()
        print("there are " + str(len(score_array)) + " reads in group " + key)

        new_score_nested_list = []
        for each_list in score_list:
            new_score_list = [0] * 324
            max_index = each_list.index(max(each_list))  # identify the MuSIC barcode index of each read
            new_score_list[max_index] = 1
            new_score_nested_list.append(new_score_list)
        sum_barcode = np.sum(new_score_nested_list, axis=0)  # sum all reads to see the distribution
        print(key, sum_barcode)

        score_of_reads_nested_list.append(sum_barcode)
        max_barcode_index = np.argmax(sum_barcode)  # identify the inferred barcode of each pMuSIC sample group
        inferred_barcode_name = ctrl_barcode[max_barcode_index][1]
        inferred_barcode_list.append(inferred_barcode_name)

    np.save(result_path + 'Step2.inferred_barcode_list.npy', inferred_barcode_list)

    # step10: calculate the percentage of each pMuSIC in 324 reference barcode
    score_of_reads_array = np.array(score_of_reads_nested_list)
    data = pd.DataFrame(score_of_reads_array)
    data.to_excel(result_path + 'Step2.distribution_of_each_pMuSIC_for_20x324_heatmap.xlsx', index=False)

    row_sums = data.sum(axis=1)
    percentage = np.round(data.div(row_sums, axis=0), decimals=3)
    print(percentage)

    percentage.to_excel(result_path + 'Step2.pMuSIC_percentage_matrix(20x324).xlsx', index=False)



