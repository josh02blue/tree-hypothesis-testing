import os
import random
import shutil
import matplotlib.pyplot as plt
from subprocess import call
from eric_functions import weighted_distance
from trees_to_edge_lengths import compute_frechet_mean, extract_frechet_mean

# ==== CONFIGURATION ====
SPECIES1_TREES = "output/gts_dendropy_CAT_tauAB-10000_tauABC-10100_tauRoot-11100_pAB-10000_pABC-10000_pRoot-10000"
SPECIES2_TREES = "output/gts_dendropy_CAT_tauAB-10000_tauABC-10100_tauRoot-11100_pAB-20000_pABC-20000_pRoot-20000"
REAL_MEAN_1 = "output/frechet_mean_gts_dendropy_CAT_tauAB-10000_tauABC-10100_tauRoot-11100_pAB-10000_pABC-10000_pRoot-10000"
REAL_MEAN_2 = "output/frechet_mean_gts_dendropy_CAT_tauAB-10000_tauABC-10100_tauRoot-11100_pAB-20000_pABC-20000_pRoot-20000"
JAR_FILE = "SturmMean_201102.jar"
OUTPUT_FILE = "output/permuted_mean_distances.txt"
NUM_PERMUTATIONS = 100

# Temporary files used per permutation
TEMP_DIR = "output/tmp_perm/"
TREE_FILE_1 = os.path.join(TEMP_DIR, "group1.tre")
TREE_FILE_2 = os.path.join(TEMP_DIR, "group2.tre")
MEAN_FILE_1 = os.path.join(TEMP_DIR, "mean1.tre")
MEAN_FILE_2 = os.path.join(TEMP_DIR, "mean2.tre")

# ==== MAIN PERMUTATION TEST ====
if __name__ == "__main__":
    os.makedirs(TEMP_DIR, exist_ok=True)

    # Load gene trees
    with open(SPECIES1_TREES, 'r') as f:
        trees1 = [line.strip() for line in f if line.strip()]
    with open(SPECIES2_TREES, 'r') as f:
        trees2 = [line.strip() for line in f if line.strip()]

    all_trees = trees1 + trees2
    labels = [0] * len(trees1) + [1] * len(trees2)

    # Load precomputed Frechet means
    mean1 = extract_frechet_mean(REAL_MEAN_1)
    mean2 = extract_frechet_mean(REAL_MEAN_2)
    observed_distance = weighted_distance(mean1, mean2)
    print(f"Test Statistic: {observed_distance:.6f}")

    distances = []

    for i in range(NUM_PERMUTATIONS):
        print("*******************************************************")
        print(f"Permutation {i}/{NUM_PERMUTATIONS}")
        print("*******************************************************")

        # Shuffle trees and relabel
        shuffled = list(zip(all_trees, labels))
        random.shuffle(shuffled)
        group1 = [tree for tree, _ in shuffled[:len(trees1)]]
        group2 = [tree for tree, _ in shuffled[len(trees1):]]

        # Write trees to files
        with open(TREE_FILE_1, 'w') as f:
            f.write("\n".join(group1) + "\n")
        with open(TREE_FILE_2, 'w') as f:
            f.write("\n".join(group2) + "\n")

        # Compute Fréchet mean of each group
        compute_frechet_mean(TREE_FILE_1, MEAN_FILE_1)
        compute_frechet_mean(TREE_FILE_2, MEAN_FILE_2)

        perm_mean1 = extract_frechet_mean(MEAN_FILE_1)
        perm_mean2 = extract_frechet_mean(MEAN_FILE_2)

        dist = weighted_distance(perm_mean1, perm_mean2)
        distances.append(dist)

    # Save all permuted distances
    with open(OUTPUT_FILE, 'w') as f:
        for d in distances:
            f.write(f"{d}\n")

    # Clean up temp folder
    shutil.rmtree(TEMP_DIR)
    print(f"Saved {NUM_PERMUTATIONS} distances to {OUTPUT_FILE}")

    # Compute and show p-value
    p_value = sum(d >= observed_distance for d in distances) / NUM_PERMUTATIONS
    print(f"P-value: {p_value:.6f}")

    # Visualization
    plt.hist(distances, bins=30, color="skyblue", edgecolor="black")
    plt.axvline(observed_distance, color="red", linestyle="--", label=f"Observed = {observed_distance:.2f}")
    plt.title("Permutation Test: Distance Between Frechet Means")
    plt.xlabel("Distance")
    plt.ylabel("Frequency")
    plt.legend()
    plt.tight_layout()
    plt.show()
