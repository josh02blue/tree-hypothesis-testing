import os, random,jpype
import matplotlib.pyplot as plt

from pathlib import Path
from eric_functions import weighted_distance
from frechet_jpype import frechet_mean_jpype # JPype Fréchet mean (SturmMean.jar)
from Pipeline.BHV_JPYPE.bhv_jpype import compute_bhv_distance   # JPype BHV distance (gtp.jar)

# ==============================
# CONFIG
# ==============================
USE_SEPARATE_FILES = False  # True: read from two files; False: split a single file
SPECIES1_TREES = "output/gts_dendropy_CAT_tauAB-10000_tauABC-10500_tauRoot-11000_pAB-10000_pABC-10000_pRoot-10000"
SPECIES2_TREES = "output/gts_dendropy_CAT_tauAB-10000_tauABC-10100_tauRoot-11100_pAB-10000_pABC-10000_pRoot-10000"

SINGLE_TREE_FILE = "rootedtree_8.txt"  # first 500 are run 1, next 500 are run 2
SPLIT_INDEX = 500

# "BHV" uses gtp.jar via JPype; "WRF" uses weighted RF function
DISTANCE_MODE = "BHV"

NUM_PERMUTATIONS = 20
OUTPUT_FILE = "output/permuted_mean_distances.txt"

# Fréchet mean controls
MEAN_NUM_ITER = 8000
MEAN_CAUCHY_LEN = 5
MEAN_EPSILON = 0.01

STURM = str(Path("SturmMean_201102.jar").resolve()) 
GTP = str(Path("gtp.jar").resolve())

# ==============================
# Distance Function
# ==============================
def distance_between(a, b):
    if DISTANCE_MODE.upper() == "BHV":
        # a and b are Java PhyloTree objects here
        return compute_bhv_distance(a, b, rooted=True)
    # for WRF we need strings:
    a_str, b_str = str(a.toString()), str(b.toString())
    if not a_str.endswith(";"): a_str += ";"
    if not b_str.endswith(";"): b_str += ";"
    return float(weighted_distance(a_str, b_str))

# ==============================
# MAIN
# ==============================
if __name__ == "__main__":
    os.makedirs("output", exist_ok=True)

    # --- start the JVM once with BOTH jars (SturmMean.jar & gtp.jar) ---
    if not jpype.isJVMStarted():
        jpype.startJVM(classpath=[STURM, GTP])

    # --- load the two groups of trees ---
    if USE_SEPARATE_FILES:
        with open(SPECIES1_TREES) as f:
            trees1 = [ln.strip() for ln in f if ln.strip()]
        with open(SPECIES2_TREES) as f:
            trees2 = [ln.strip() for ln in f if ln.strip()]
    else:
        with open(SINGLE_TREE_FILE) as f:
            all_lines = [ln.strip() for ln in f if ln.strip()]
        trees1 = all_lines[:SPLIT_INDEX]
        trees2 = all_lines[SPLIT_INDEX:SPLIT_INDEX*2]

    # --- compute observed Fréchet means and their distance (computing our test statistic) ---
    mean1 = frechet_mean_jpype(trees1, num_iter=MEAN_NUM_ITER, cauchy_len=MEAN_CAUCHY_LEN,epsilon=MEAN_EPSILON, return_java=True)
    mean2 = frechet_mean_jpype(trees2, num_iter=MEAN_NUM_ITER, cauchy_len=MEAN_CAUCHY_LEN,epsilon=MEAN_EPSILON, return_java=True)
    observed_distance = distance_between(mean1, mean2)
    print(f"\nDistance mode: {DISTANCE_MODE}")
    print(f"Observed distance: {observed_distance:.6f}\n")
    # ==============================
    # permutation setup
    # ==============================
    # --- Assigning labels ---
    all_trees = trees1 + trees2
    labels = [0] * len(trees1) + [1] * len(trees2)
    n1 = len(trees1)

    # --- run permutations: recompute means for shuffled groups, then distance ---
    distances = []
    for i in range(1, NUM_PERMUTATIONS + 1):
        print(f"--- Permutation {i}/{NUM_PERMUTATIONS} ---")

        shuffled = list(zip(all_trees, labels))
        random.shuffle(shuffled)
        group1 = [t for t, _ in shuffled[:n1]]
        group2 = [t for t, _ in shuffled[n1:]]

        perm_mean1 = frechet_mean_jpype(group1, num_iter=MEAN_NUM_ITER, cauchy_len=MEAN_CAUCHY_LEN,epsilon=MEAN_EPSILON, return_java=True)
        perm_mean2 = frechet_mean_jpype(group2, num_iter=MEAN_NUM_ITER, cauchy_len=MEAN_CAUCHY_LEN,epsilon=MEAN_EPSILON, return_java=True)

        dist = distance_between(perm_mean1, perm_mean2)
        distances.append(dist)
        print(f"  perm distance: {dist:.6f}")

    # --- save distances so they can be pooled across runs/people ---
    with open(OUTPUT_FILE, "w") as f:
        for d in distances:
            f.write(f"{d}\n")
    print(f"\nSaved {NUM_PERMUTATIONS} distances → {OUTPUT_FILE}")

    # --- compute simple permutation p-value ---
    p_value = sum(d >= observed_distance for d in distances) / max(1, NUM_PERMUTATIONS)
    print(f"Permutation p-value: {p_value:.6f}")

    # --- plot histogram with observed distance marker ---
    plt.hist(distances, bins=30, color="skyblue", edgecolor="black")
    plt.axvline(observed_distance, color="red", linestyle="--", label=f"Observed = {observed_distance:.2f}")
    plt.title("Permutation Test: Distance Between Frechet Means")
    plt.xlabel("Distance")
    plt.ylabel("Frequency")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # --- optional: shut down JVM when done ---
    if jpype.isJVMStarted():
        jpype.shutdownJVM()
