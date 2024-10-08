import numpy as np
from Bio import Align
from Bio.Seq import Seq
import pandas as pd

project_root = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join(project_root, 'output/')
# step6: create the control library

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

fp = []
for i in range(1, 19):
    name = "fp_"+str(i)
    fp.append(eval(name))  # change the str "fp_1" to variable fp_1
    fp[i-1] = fp[i-1].upper()
print("the length of fp1 is " + str(len(fp[0])))
print("the sequence of fp1 is:")
print(fp[0])
print("the length of list fp is " + str(len(fp)))

cmv_fp = []  # list "cmv_fp" contains NheI-CMV-BamHI-Kozak-ATG-fp
for i in range(len(fp)):
    sequence = NheI_cmv_kozak.upper() + fp[i]
    cmv_fp.append(sequence)
print("the length of NheI_cmv_BamHI_kozak_fp1 is:")
print(len(cmv_fp[0]))

fp = np.array(fp)
np.save(result_path + 'Step2_fp.npy', fp)

cmv_fp = np.array(cmv_fp) # expected read without the sequence tail to MfeI site
np.save(result_path + 'Step2_cmv_fp.npy', cmv_fp)

# step7: align and score to cmv_fp to minimize the global scoring errors
final_good_read_pool = np.load(result_path + 'Step1_good_alignment_pool.npy', allow_pickle=True).item()

aligner = Align.PairwiseAligner()
aligner.mode = "global"

score_of_read = {}
for key in list(final_good_read_pool.keys()):
    fp_list = final_good_read_pool[key]
    fp_score_list = []
    for each_seq in fp_list:
        score_list = []
        for each_ctrl in cmv_fp:
            each_score = aligner.score(Seq(each_seq[:-2]), Seq(each_ctrl))
            score_list.append(each_score)
        fp_score_list.append(score_list)
    fp_score_array = np.array(fp_score_list)
    score_of_read.update({key: fp_score_array})
    print('the scoring results for ' + key + " are:")
    print('the size of fp_score_array is ', fp_score_array.size)
    print('the shape of fp_score_array is ', fp_score_array.shape)

np.save("score_of_read.npy", score_of_read, allow_pickle=True)
print("now you have score_of_read stored!")

# step8: identify each read
fp_name = ["EBFP2", "mTagBFP2", "mT-Sapphire", "mAmetrine", "mCerulean3", "LSSmOrange", "mBeRFP", "mTFP1", "EGFP",
           "CyOFP1", "mClover3", "mVenus", "mPapaya", "mOrange2", "mRuby3", "mKate2", "mCardinal", "miRFP670"]

score_of_reads_nested_list = []
inferred_fp_list = []

for key in list(score_of_read.keys()):
    score_array = score_of_read[key]
    score_list = score_array.tolist()
    print("there are " + str(len(score_array)) + " reads in group " + key)
    new_score_nested_list = []
    for each_list in score_list:
        new_score_list = [0] * 18
        max_index = each_list.index(max(each_list))
        # identify the inferred probe of each read
        new_score_list[max_index] = 1
        new_score_nested_list.append(new_score_list)
    sum_fp = np.sum(new_score_nested_list, axis=0)  # along the column
    print(key, sum_fp)

    # identify the inferred probe of each barcode file
    max_probe_index = np.argmax(sum_fp)
    inferred_fp_name = fp_name[max_probe_index]
    inferred_fp_list.append(inferred_fp_name)
    score_of_reads_nested_list.append(sum_fp)
print(inferred_fp_list)
score_of_reads_array = np.array(score_of_reads_nested_list)

np.save(result_path + 'Step2.inferred_fp_list.npy', inferred_fp_list)
np.save(result_path + 'Step2.actual_fp_list.npy', fp_name)

data = pd.DataFrame(score_of_reads_array)
data.to_excel(result_path + 'Step2.distribution_of_each_barcode_for_18x18_heatmap.xlsx', index=False)

# step9: calculate the true positive rate
row_sums = data.sum(axis=1)
percentage = np.round(data.div(row_sums, axis=0), decimals=2)

score_of_reads_array = np.sum(score_of_reads_array, axis=0)
print(score_of_reads_array)

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
print(percentage)

percentage.to_excel(result_path + 'Step2.final_percentage_for_18x18_heatmap.xlsx', index=False)
