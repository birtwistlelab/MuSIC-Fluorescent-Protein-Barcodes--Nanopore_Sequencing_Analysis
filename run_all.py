import subprocess
import argparse


def run_step(folder, step):
    script_path = f"./{folder}/{step}.py"
    print(f"Running {script_path}...") 
    subprocess.run(["python", script_path], check=True)


def main():
    parser = argparse.ArgumentParser(description="Run experiemnts.")
    parser.add_argument('--experiment', choices=['exp1', 'exp2', 'exp3', 'all'], default='all',
                        help="Specify which experiment to run ('exp1', 'exp2', 'exp3' or 'all').")
    parser.add_argument('--step', choices=['Step1', 'Step2', 'Step3', 'all'], default='all',
                        help="Specify which step to run('Step1', 'Step2', 'Step3', or 'all').")
    args = parser.parse_args()

    experiments = {
        'exp1': '1.pR-fp_pool',
        'exp2': '2.sequenced_pMuSIC_pool',
        'exp3': '3.MuSIC_barcodes_pool',
    }

    step_mapping = {
        'Step1': 'Step1.Sequence_Processing',
        'Step2': 'Step2.Alignment_and_Analysis',
        'Step3': 'step3.Heatmap_Generation',
    }

    if args.experiment == 'all':
        for folder in experiments.values():
            run_all_steps(folder, args.step, step_mapping)
    else:
        run_all_steps(experiments[args.experiment], args.step, step_mapping)


def run_all_steps(folder, step, step_mapping):
    if step == 'all':
        for step_name in ['Step1', 'Step2', 'Step3']:
            run_step(folder, f"{folder.split('.')[0]}_{step_mapping[step_name]}")
    else:
        run_step(folder, f"{folder.split('.')[0]}_{step_mapping[step]}")
    

if __name__ == "__main__":
    main()

