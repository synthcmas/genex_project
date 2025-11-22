#!/usr/bin/env python3
import argparse
import os
import math
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROTEIN_DEGRADATION_RATE = 0.0005

# -------------------- Utilities --------------------
def round_up_to_nearest_10(x):
    return int(math.ceil(x / 10.0)) * 10

def get_num_from_datafile(path, default):
    if not os.path.exists(path):
        return default
    try:
        values = []
        with open(path, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    # second column is the count/value used earlier
                    values.append(float(parts[1]))
        if not values:
            return default
        max_val = max(values)
        return round_up_to_nearest_10(max_val) + 20
    except Exception:
        return default

def scale_parameters(params):
    """
    Input: unscaled [nu0, nu1, d0, d1, k0, k1]
    Output: scaled [aa, bb, gamma, kappa0, kappa1]
    """
    nu0, nu1, d0, d1, k0, k1 = params
    aa = nu0 / d1
    bb = nu1 / d0
    gamma = d0 / d1
    kappa0 = k0 / d1
    kappa1 = k1 / d1
    return [aa, bb, gamma, kappa0, kappa1]

# -------------------- Core Logic --------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Python wrapper for genex")
    parser.add_argument("-p", nargs=5, type=float, required=True,
                        help="Five parameters (scaled or unscaled)")
    parser.add_argument("-s", action="store_true",
                        help="If set, parameters are already scaled")
    return parser.parse_args()

def prepare_parameters(params, scaled):
    if scaled:
        return params  # already scaled
    else:
        nu0, nu1, d0, k0, k1 = params
        d1 = PROTEIN_DEGRADATION_RATE
        return scale_parameters([nu0, nu1, d0, d1, k0, k1])

def create_parameters_ini(aa, bb, gamma, kappa0, kappa1, num_mrna, num_protein):
    return f"""# Number of mRNA appearing in the master equation
Num_mRNA={num_mrna}

# Number of Proteins appearing in the master equation
Num_Protein={num_protein}


# The probability of transcription per unit time
#Transcription_Probability=0.002
           
aa = {aa}
          
# This is nu0/d1
# The probability of translation per unit time
#Translation_Probability=0.05
bb = {bb}
# This is nu1/d0
# The probability of degradation of mRNA per unit time
#mRNA_Degradation_Probability=0.005
# The probability of degradation of proteins per unit time
#Protein_Degradation_Probability=0.0005
gamma = {gamma}
# This is d0/d1
# The probability of transition from inactive promoter to active promoter
#Inactive_to_Active_Transition_probability=0.0003
kappa0 = {kappa0}
# k0/d1
# The probability of transition from active promoter to inactive promoter
#Active_to_Inactive_Transition_probability=0.0001
kappa1 = {kappa1}

#k1/d1
# Number of snapshots in time needed to track the evolution of the probabilities
SnapShots=2


# Initial time
Initial_Time=1
# Final time 
Final_Time=1000
"""

def write_parameters_ini(content, filename="parameters.ini"):
    with open(filename, "w") as f:
        f.write(content)
    print(f"✅ {filename} written successfully.")

def run_genex(filename="parameters.ini"):
    try:
        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        fortran_dir = os.path.join(project_root, "fortran")
        abs_filename = os.path.join(project_root, filename)

        print(f"🛠 Compiling genex in {fortran_dir} ...")
        subprocess.run(["make", "clean"], cwd=fortran_dir, check=True)
        subprocess.run(["make"], cwd=fortran_dir, check=True)
        print("✅ Compilation successful.")

        exe_path = os.path.join(fortran_dir, "genex")

        # 👇 Run Fortran from project root where parameters.ini exists
        subprocess.run([exe_path, filename], cwd=project_root, check=True)
        print("✅ genex executed successfully.")

    except FileNotFoundError as e:
        print(f"⚠️ Command not found: {e}")
    except subprocess.CalledProcessError as e:
        print(f"⚠️ Command failed: {e}")


# -------------------- Data & Plotting --------------------
def load_data(mrna_file, protein_file, mRNA_data_file, protein_data_file):
    df_mRNA = pd.read_csv(mrna_file, delim_whitespace=True)
    df_protein = pd.read_csv(protein_file, delim_whitespace=True)

    mRNA_bin_start, mRNA_bin_stop, mRNA_prob, mRNA_err = np.loadtxt(
        mRNA_data_file, usecols=(0,1,2,3), unpack=True
    )
    protein_bin_start, protein_bin_stop, protein_prob, protein_err = np.loadtxt(
        protein_data_file, usecols=(0,1,2,3), unpack=True
    )
    return df_mRNA, df_protein, (mRNA_bin_start, mRNA_bin_stop, mRNA_prob, mRNA_err), (protein_bin_start, protein_bin_stop, protein_prob, protein_err)

def compute_binned_theory(df, bin_start, bin_stop):
    theory = np.zeros_like(bin_start)
    for i in range(len(bin_start)):
        start, stop = int(bin_start[i]), int(bin_stop[i])
        theory[i] = sum(df.Probability[start:stop]) / sum(df.Probability)
    return theory / sum(theory)

def plot_distribution(bin_mean, theory, data_prob, data_err, xlabel, filename, param_text):
    plt.style.use('seaborn-darkgrid')
    fig, ax = plt.subplots(figsize=(7,7), dpi=80, facecolor='#f7f7f7')
    ax.set_facecolor('white')

    # Plot theory and data
    ax.plot(bin_mean, theory, lw=5.0, color='#2a9d8f', label='3-Stage')
    ax.errorbar(bin_mean, data_prob, yerr=data_err, fmt='o', linestyle="", c='#e76f51', label='Data', lw=3)

    # Grid, labels, ticks
    ax.grid(True, linestyle='--', linewidth=0.6, alpha=0.7)
    ax.set_xlabel(xlabel, fontsize=22, fontweight='bold')
    ax.set_ylabel('Probability', fontsize=22, fontweight='bold')
    ax.tick_params(axis='both', which='both', color='black', labelcolor='black', direction='out', length=6, width=1.2)
    ax.legend(fontsize=24, frameon=True, facecolor='white', edgecolor='black', framealpha=1, shadow=False, loc='best')

    # Parameter box
    ax.text(
        0.9, 0.6, param_text,
        transform=ax.transAxes, fontsize=18,
        verticalalignment='top', horizontalalignment='right',
        color='#222222',
        bbox=dict(boxstyle="round,pad=0.7", facecolor=(1,1,1,0.85), edgecolor='black', linewidth=1.5)
    )

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1.2)

    plt.tight_layout()
    plt.savefig(filename)
    plt.show()
    print(f"✅ Saved figure to {filename}")

# -------------------- Run Wrapper --------------------
def run_pipeline(params, scaled=False):
    """
    Run full pipeline (INI creation → genex → plotting).
    params: list/tuple of 5 floats (scaled or unscaled depending on scaled flag)
    scaled: bool, True if params are already scaled
    """
    aa, bb, gamma, kappa0, kappa1 = prepare_parameters(params, scaled)
    num_mrna = get_num_from_datafile("./data/mRNAdata.txt", 150)
    num_protein = get_num_from_datafile("./data/proteindata.txt", 50)

    ini_content = create_parameters_ini(aa, bb, gamma, kappa0, kappa1, num_mrna, num_protein)
    write_parameters_ini(ini_content)
    run_genex()

    df_mRNA, df_protein, mRNA_bins, protein_bins = load_data(
        "plot_data/Genex3_OUT_mRNA_timeslice_2.txt",
        "plot_data/Genex3_OUT_Protein_timeslice_2.txt",
        "data/mRNAdata.txt",
        "data/proteindata.txt"
    )

    mRNA_bin_mean = (mRNA_bins[0] + mRNA_bins[1]) / 2
    protein_bin_mean = (protein_bins[0] + protein_bins[1]) / 2
    mRNA_theory = compute_binned_theory(df_mRNA, mRNA_bins[0], mRNA_bins[1])
    protein_theory = compute_binned_theory(df_protein, protein_bins[0], protein_bins[1])

    param_text = (
        r"$\mathtt{{a\ =\ {:>5.1f}}}$" "\n"
        r"$\mathtt{{b\ =\ {:>5.2f}}}$" "\n"
        r"$\mathtt{{\gamma\ =\ {:>5.2f}}}$" "\n"
        r"$\mathtt{{\kappa_0\ =\ {:>5.2f}}}$" "\n"
        r"$\mathtt{{\kappa_1\ =\ {:>5.2f}}}$"
    ).format(aa, bb, gamma, kappa0, kappa1)

    plot_distribution(mRNA_bin_mean, mRNA_theory, mRNA_bins[2], mRNA_bins[3], 'mRNA', 'figs/mRNA_Data_Compare.png', param_text)
    plot_distribution(protein_bin_mean, protein_theory, protein_bins[2], protein_bins[3], 'Protein', 'figs/Protein_Data_Compare.png', param_text)

# -------------------- CLI Entry Point --------------------
def main():
    args = parse_args()
    run_pipeline(args.p, args.s)

if __name__ == "__main__":
    main()
