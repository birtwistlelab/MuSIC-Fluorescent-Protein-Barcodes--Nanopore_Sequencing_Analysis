## **Installation**
For specific Python packages, a requirements file has been provided at 'requirements.txt'.
Assuming Pythonis installed, run the following command to install all necessary dependencies:
```
pip install -r requirements.txt
```

## **Usage**
All result files were available upon request due to size limitations.

To use the project, you can use the run_all.py file in the terminal for both top-level direction for complete functionality, or you can test each experiment individually for specific tasks.
All experiments and their steps are organized as shown in the structure tree below: 


![structure tree](structure_tree.png)


## **Running the Full Operation**
```
python run_all.py --experiment all --step all
```

### **Running Individual Experiment**
```
python run_all.py --experiment exp1  --step all
python run_all.py --experiment exp2  --step all
python run_all.py --experiment exp3  --step all
```

### **Running Individual Scripts**
```
python run_all.py --experiment exp1  --step Step1
python run_all.py --experiment exp1  --step Step2
python run_all.py --experiment exp1  --step Step3
```

To run other experiments, simply replace exp1 with the desired experiment name, such as exp2 or exp3. 
To view Step 2 or Step 3 of one experiment, you should first run the previous step(s); otherwise, you won't be able to directly access the results of the prior step locally.

## **Warning**
The dataset used in this project is substantial, and certain steps may require extended processing time. 
For instance, on a 2 GHz Quad-Core Intel Core i5 with 16 GB 3733 MHz LPDDR4X:
Step 2 of Experiment 2 takes more than 12 hours to complete.
Step 2 of Experiment 3 takes approximately 36 hours to complete.
It is **strongly recommended** to run this code on a **High-Performance Computing (HPC) cluster**.
