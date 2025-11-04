import os
import csv
import random
from pathlib import Path
import jpype
import matplotlib
matplotlib.use("Agg")  
import matplotlib.pyplot as plt


from frechet_jpype import frechet_mean_jpype   # returns Java PhyloTree if return_java=True
from Pipeline.BHV_JPYPE.bhv_jpype import compute_bhv_distance     # can take Java PhyloTree objects directly

# =========================
# CONFIG (edit as needed)
# =========================
INPUT_DIR        = "rootedtrees"   # folder containing rootedtree_#.txt files
FILE_PREFIX      = "rootedtree_"
START_INDEX      = 1                 # inclusive
END_INDEX        = 304               # inclusive 
SPLIT_INDEX      = 500               # first 500 vs next 500 trees per file

# Distance mode: "BHV" (via gtp.jar) or "WRF" (your weighted RF)
DISTANCE_MODE    = "BHV"

# Freechet mean parameters 
MEAN_NUM_ITER    = 8000
MEAN_CAUCHY_LEN  = 5
MEAN_EPSILON     = 0.01

# Permutation settings
NUM_PERMUTATIONS = 50                # adjust if you want more/less
RANDOM_SEED      = 123               # make runs reproducible

# Outputs
OUTPUT_DIR       = "batch_output"
CSV_PATH         = os.path.join(OUTPUT_DIR, "new_permutation_summary.csv")
SAVE_PLOTS       = True              # set False if you don’t want PNG histograms
PLOTS_DIR        = os.path.join(OUTPUT_DIR, "plots")

# JARs 
STURM_JAR        = str(Path("SturmMean_201102.jar").resolve())
GTP_JAR          = str(Path("gtp.jar").resolve())

# =========================
# Helpers
# =========================
def ensure_jvm():
    """Start JVM once with both JARs on the classpath."""
    if not jpype.isJVMStarted():
        jpype.startJVM(classpath=[STURM_JAR, GTP_JAR])

def distance_between(a, b):
    """
    Distance between two mean trees
    If BHV: a and b are Java PhyloTree objects (fast path)
    If WRF: convert to Newick strings for weighted RF
    """
    if DISTANCE_MODE.upper() == "BHV":
        return float(compute_bhv_distance(a, b, rooted=True))
    from eric_functions import weighted_distance
    a_str, b_str = str(a.toString()), str(b.toString())
    if not a_str.endswith(";"): a_str += ";"
    if not b_str.endswith(";"): b_str += ";"
    return float(weighted_distance(a_str, b_str))

def run_permutation_for_file(tree_path):
    """
    Run the whole pipeline for one locus file
    Returns a dict with:
      filename, n_perm, null_min, null_max, null_mean, test_stat, p_value
    """
    # load and split trees
    with open(tree_path) as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    trees1 = lines[:SPLIT_INDEX]
    trees2 = lines[SPLIT_INDEX:SPLIT_INDEX*2]

    # compute observed means 
    mean1 = frechet_mean_jpype(
        trees1, num_iter=MEAN_NUM_ITER, cauchy_len=MEAN_CAUCHY_LEN,
        epsilon=MEAN_EPSILON, return_java=True)
    mean2 = frechet_mean_jpype(
        trees2, num_iter=MEAN_NUM_ITER, cauchy_len=MEAN_CAUCHY_LEN,
        epsilon=MEAN_EPSILON, return_java=True)

    test_stat = distance_between(mean1, mean2)

    # permutations
    all_trees = trees1 + trees2
    labels = [0]*len(trees1) + [1]*len(trees2)
    n1 = len(trees1)
    null_distances = []

    for i in range(NUM_PERMUTATIONS):
        shuffled = list(zip(all_trees, labels))
        random.shuffle(shuffled)
        group1 = [t for t,_ in shuffled[:n1]]
        group2 = [t for t,_ in shuffled[n1:]]

        pmean1 = frechet_mean_jpype(
            group1, num_iter=MEAN_NUM_ITER, cauchy_len=MEAN_CAUCHY_LEN,
            epsilon=MEAN_EPSILON, return_java=True)
        pmean2 = frechet_mean_jpype(
            group2, num_iter=MEAN_NUM_ITER, cauchy_len=MEAN_CAUCHY_LEN,
            epsilon=MEAN_EPSILON, return_java=True)

        null_distances.append(distance_between(pmean1, pmean2))

    # stats
    null_min  = min(null_distances) if null_distances else float("nan")
    null_max  = max(null_distances) if null_distances else float("nan")
    null_mean = (sum(null_distances)/len(null_distances)) if null_distances else float("nan")
    p_value   = (sum(d >= test_stat for d in null_distances) / max(1, len(null_distances)))

    # optional plot
    if SAVE_PLOTS:
        os.makedirs(PLOTS_DIR, exist_ok=True)
        base = Path(tree_path).stem
        out_png = os.path.join(PLOTS_DIR, f"{base}_null_hist.png")
        plt.figure()
        plt.hist(null_distances, bins=30)
        plt.axvline(test_stat, linestyle="--", label=f"Observed = {test_stat:.3f}")
        plt.title(f"{base}: Permutation Null ({DISTANCE_MODE})")
        plt.xlabel("Distance between perm means"); plt.ylabel("Frequency")
        plt.legend(); plt.tight_layout()
        plt.savefig(out_png, dpi=150)
        plt.close()

    return {
        "filename": Path(tree_path).name,
        "n_perm": len(null_distances),
        "null_min": null_min,
        "null_max": null_max,
        "null_mean": null_mean,
        "test_stat": test_stat,
        "p_value": p_value,
    }

# =========================
# Main
# =========================
if __name__ == "__main__":
    random.seed(RANDOM_SEED)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    ensure_jvm()

    # Prepare CSV (append if exists, else write header)
    write_header = not os.path.exists(CSV_PATH)
    with open(CSV_PATH, "a", newline="") as csvfile:
        writer = csv.DictWriter(
            csvfile,
            fieldnames=["filename", "n_perm", "null_min", "null_max", "null_mean", "test_stat", "p_value"]
        )
        if write_header:
            writer.writeheader()

        # Iterate files
        for i in range(START_INDEX, END_INDEX + 1):
            tree_path = os.path.join(INPUT_DIR, f"{FILE_PREFIX}{i}.txt")
            if not os.path.exists(tree_path):
                print(f"[skip] {tree_path} not found")
                continue

            print(f"[run ] {Path(tree_path).name}")
            try:
                result = run_permutation_for_file(tree_path)
                writer.writerow(result)
                csvfile.flush()
                print(f"[done] {result['filename']} | stat={result['test_stat']:.4f} p={result['p_value']:.4f}")
            except Exception as e:
                print(f"[err ] {Path(tree_path).name}: {e}")

    # Optional: shut down JVM
    if jpype.isJVMStarted():
        jpype.shutdownJVM()