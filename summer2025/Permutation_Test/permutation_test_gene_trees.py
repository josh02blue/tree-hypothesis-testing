import os
import random
import shutil
import matplotlib.pyplot as plt
from eric_functions import weighted_distance
from trees_to_edge_lengths import compute_frechet_mean, extract_frechet_mean
from Pipeline.BHV_JPYPE.bhv_jpype import compute_bhv_distance

# ==== CONFIGURATION ====
# === OPTION A: Use two separate gene tree files ===
USE_SEPARATE_FILES = False  # Set to True to use the two SPECIES files instead
USE_BHV_DISTANCE = False     # Set to True to use BHV instead of RF distance

# If using separate gene tree files
SPECIES1_TREES = "output/gts_dendropy_CAT_tauAB-10000_tauABC-10100_tauRoot-11100_pAB-10000_pABC-10000_pRoot-10000"
SPECIES2_TREES = "output/gts_dendropy_CAT_tauAB-10000_tauABC-10100_tauRoot-11100_pAB-20000_pABC-20000_pRoot-20000"
REAL_MEAN_1 = "output/frechet_mean_gts_dendropy_CAT_tauAB-10000_tauABC-10100_tauRoot-11100_pAB-10000_pABC-10000_pRoot-10000"
REAL_MEAN_2 = "output/frechet_mean_gts_dendropy_CAT_tauAB-10000_tauABC-10100_tauRoot-11100_pAB-20000_pABC-20000_pRoot-20000"

# === OPTION B: Use one .txt file and split it ===
SINGLE_TREE_FILE = "rootedtree_7.txt"
SPLIT_INDEX = 500
 
OUTPUT_FILE = "output/permuted_mean_distances.txt"
NUM_PERMUTATIONS = 50

TEMP_DIR = "output/tmp_perm/"
TREE_FILE_1 = os.path.join(TEMP_DIR, "group1.tre")
TREE_FILE_2 = os.path.join(TEMP_DIR, "group2.tre")
MEAN_FILE_1 = os.path.join(TEMP_DIR, "mean1.tre")
MEAN_FILE_2 = os.path.join(TEMP_DIR, "mean2.tre")

if __name__ == "__main__":
    os.makedirs(TEMP_DIR, exist_ok=True)

    if USE_SEPARATE_FILES:
        with open(SPECIES1_TREES, 'r') as f:
            trees1 = [line.strip() for line in f if line.strip()]
        with open(SPECIES2_TREES, 'r') as f:
            trees2 = [line.strip() for line in f if line.strip()]
    else:
        with open(SINGLE_TREE_FILE, 'r') as f:
            all_lines = [line.strip() for line in f if line.strip()]
        trees1 = all_lines[:SPLIT_INDEX]
        trees2 = all_lines[SPLIT_INDEX:SPLIT_INDEX*2] 

    all_trees = trees1 + trees2
    labels = [0] * len(trees1) + [1] * len(trees2)

    if USE_SEPARATE_FILES:
        mean1 = extract_frechet_mean(REAL_MEAN_1)
        mean2 = extract_frechet_mean(REAL_MEAN_2)
    else:
        with open("output/run1_split.tre", "w") as f:
            f.write("\n".join(trees1) + "\n")
        with open("output/run2_split.tre", "w") as f:
            f.write("\n".join(trees2) + "\n")

        compute_frechet_mean("output/run1_split.tre", "output/mean_run1.tre")
        compute_frechet_mean("output/run2_split.tre", "output/mean_run2.tre")

        mean1 = extract_frechet_mean("output/mean_run1.tre")
        mean2 = extract_frechet_mean("output/mean_run2.tre")

    # Compute observed distance
    if USE_BHV_DISTANCE:
        observed_distance = compute_bhv_distance(mean1, mean2)
    else:
        observed_distance = weighted_distance(mean1, mean2)

    print(f"Test Statistic: {observed_distance:.6f}")
    distances = []

    for i in range(NUM_PERMUTATIONS):
        print("*******************************************************")
        print(f"Permutation {i}/{NUM_PERMUTATIONS}")
        print("*******************************************************")

        shuffled = list(zip(all_trees, labels))
        random.shuffle(shuffled)
        group1 = [tree for tree, _ in shuffled[:len(trees1)]]
        group2 = [tree for tree, _ in shuffled[len(trees1):]]

        with open(TREE_FILE_1, 'w') as f:
            f.write("\n".join(group1) + "\n")
        with open(TREE_FILE_2, 'w') as f:
            f.write("\n".join(group2) + "\n")

        compute_frechet_mean(TREE_FILE_1, MEAN_FILE_1)
        compute_frechet_mean(TREE_FILE_2, MEAN_FILE_2)

        perm_mean1 = extract_frechet_mean(MEAN_FILE_1)
        perm_mean2 = extract_frechet_mean(MEAN_FILE_2)

        if USE_BHV_DISTANCE:
            dist = compute_bhv_distance(perm_mean1, perm_mean2)
        else:
            dist = weighted_distance(perm_mean1, perm_mean2)

        distances.append(dist)

    with open(OUTPUT_FILE, 'w') as f:
        for d in distances:
            f.write(f"{d}\n")

    shutil.rmtree(TEMP_DIR)
    print(f"Saved {NUM_PERMUTATIONS} distances to {OUTPUT_FILE}")

    p_value = sum(d >= observed_distance for d in distances) / NUM_PERMUTATIONS
    print(f"P-value: {p_value:.6f}")

    plt.hist(distances, bins=30, color="skyblue", edgecolor="black")
    plt.axvline(observed_distance, color="red", linestyle="--", label=f"Observed = {observed_distance:.2f}")
    plt.title("Permutation Test: Distance Between Frechet Means")
    plt.xlabel("Distance")
    plt.ylabel("Frequency")
    plt.legend()
    plt.tight_layout()
    plt.show()

