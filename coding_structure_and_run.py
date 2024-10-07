import tkinter as tk
from tkinter import messagebox
import subprocess
import matplotlib.pyplot as plt
import numpy as np

experiments = {
    'Experiment 1': './1.pR-fp_pool/1_Step1.Sequence_Processing.py',
    'Experiment 2': './2.sequenced_pMusIC_pool/2_Step1.Sequence_Processing.py',
    'Experiment 3': './3.MuSIC_barcodes_pool/3_Step1.Sequence_Processing.py',
}


def run_experiment(script):
  try:
    subprocess.run(['python', script], check=True)
    messagebox.showinfo('Success', f'Successfully ran: {script}')
  except subprocess.CalledProcessError as e:
    messagebox.showerror('Error, f'Error running {script}: {e}')


def create_tree_structure():
  fig, ax = plt.subplots(figsize=(5, 3))
  y_pos = np.arange(len(experiments))
  ax.barh(y_pos, [1] * len(experiments), align='center', alpha=0.5)
  ax.set_yticks(y_pos)
  ax.set_yticklabels(list(experiments.keys()))
  ax.invert_yaxis()
  ax.set_xlabel('Experiments')
  plt.show()


root = tk.Tk()
root.title("coding_structure_and_run")

for exp, script in experiments.items():
    btn = tk.Button(root, text=f"Run {exp}", command=lambda s=script: run_experiment(s))
    btn.pack(pady=10)


tree_button = tk.Button(root, text="Show Experiment Structure", command=create_tree_structure)
tree_button.pack(pady=20)

root.mainloop()
  
    
