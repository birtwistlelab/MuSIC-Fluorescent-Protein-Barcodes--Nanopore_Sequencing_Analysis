Installation
For specific Python packages, a requirements file has been provided at 'requirements.txt'.
Assuming Pythonis installed, run the following command to install all necessary dependencies:
pip install -r requirements.txt

Usage
All original files were uploaded via Git LFS because of the size limitation, so you can view and work with these files directly in Github Codespaces.
If you want to view these original files locally, please make sure you have installed and configuresd Git LFS.

To use the project, you can use the run_all.py file in the terminal for both top-level direction for complete functionality, or you can test each experiment individually for specific tasks.
All experiments and their steps are organized as shown in the structure tree. 
![structure tree](structure_tree.png)

Running the Full Operation
python run_all.py --experiment all --step all

Running Individual Experiment
python run_all.py --experiment exp1  --step all
python run_all.py --experiment exp2  --step all
python run_all.py --experiment exp3  --step all

Running Individual Scripts
python run_all.py --experiment exp1  --step Step1
python run_all.py --experiment exp1  --step Step2
python run_all.py --experiment exp1  --step Step3

To run other experiments, simply replace exp1 with the desired experiment name, such as exp2 or exp3.
