import random
import numpy as np
import matplotlib.pyplot as plt
import jpype
from pathlib import Path
from Pipeline.BHV_JPYPE.bhv_jpype import compute_bhv_distance
from frechet_jpype import frechet_mean_jpype

# ====== JAR setup ======
# Resolve absolute paths for Java JARs that handle BHV distances & Fréchet means
GTP_JAR   = str(Path("gtp.jar").resolve())
STURM_JAR = str(Path("SturmMean_201102.jar").resolve())

# Start JVM with the required JARs if it hasn’t been started yet
if not jpype.isJVMStarted():
    jpype.startJVM(classpath=[GTP_JAR, STURM_JAR])

# ====== USER CONFIGURATION ======
PATH        = "./rootedtrees/rootedtree_136.txt"
MAX_TREES   = 100
ITERS       = 75
EPSILON     = 0.0001
SEED        = 21
MEAN_ITERS  = 8000
MAX_K       = 8

# Optionally track cluster reassignments per iteration for diagnostics
TRACK_CHANGES = True
K_FOR_CHANGE_PLOT = 2     # which k to track when TRACK_CHANGES = True
SEED_FOR_CHANGE_PLOT = 275
# ========================

# Set random seeds for reproducibility
random.seed(SEED); np.random.seed(SEED)



# ====== Helpers ======
def ensure_semicolon(s: str) -> str:
    #Ensure the Newick string ends with a semicolon.

    s = s.strip()
    return s if s.endswith(";") else s + ";"
    

def clean_newick(s: str) -> str:
    """
    
    Extract valid Newick subtree between the first '(' and last ')',
    then ensure it ends with ';'.
    This is helpful when the Java mean output contains metadata or formatting noise.
    """
    s = str(s)
    i = s.find("(")
    j = s.rfind(")")
    if i == -1 or j == -1 or j <= i:
        raise ValueError(f"No valid Newick parentheses found in mean output; {s}")
    return ensure_semicolon(s[i:j+1])

def bhv_safe(a: str, b: str) -> float:
    """
    Compute BHV distance between two trees with fallbacks:
    - try raw strings
    - if fails, clean one or both strings
    """
    try:
        return compute_bhv_distance(a, b, rooted=True)
    except Exception:
        try:
            return compute_bhv_distance(clean_newick(a), ensure_semicolon(b), rooted=True)
        except Exception:
            # last try: clean both
            return compute_bhv_distance(clean_newick(a), clean_newick(b), rooted=True)

# ====== Load trees from file ======

with open(PATH) as f:
    trees = [ensure_semicolon(ln.strip()) for ln in f if ln.strip()]

# Sanity checks
if len(trees) < 2:
    raise SystemExit("Need at least 2 trees.")
if len(trees) > MAX_TREES:
    trees = trees[:MAX_TREES]
n = len(trees)
print(f"Using {n} trees")

# ====== k-means in BHV tree space with early stop on assignments ======
def kmeans_frechet(
    trees, k, iters=200, epsilon=0.01, seed=None, mean_iters=2000, track_changes=False
):
    """
    Perform k-means clustering in BHV tree space using Fréchet means as centroids.

    Args:
        trees: list of Newick strings
        k: number of clusters
        iters: max k-means iterations
        epsilon: tolerance for Fréchet mean convergence
        seed: RNG seed
        mean_iters: iterations for Fréchet mean computation
        track_changes: if True, record how many trees switch clusters each iteration

    Returns:
        inertia: final within-cluster sum of squares (BHV distance squared)
        iters_used: number of iterations actually run
        converged: True if labels or centers stopped changing before max iters
        changes_per_iter: list of how many trees changed clusters each iteration
    """
    # Reseed RNG if provided
    if seed is not None:
        random.seed(seed); np.random.seed(seed)
    # Randomly pick initial cluster centers
    centers = random.sample(trees, k)
    prev_labels = None
    converged = False
    changes_per_iter = []  # filled only if track_changes=True

    for it in range(1, iters + 1):
        # --- Assignment step ---
        clusters = [[] for _ in range(k)]
        lab_list = []
        for t in trees:
            # compute BHV distance to each center, assign to closest
            dists = [bhv_safe(c, t) for c in centers]
            j = int(np.argmin(dists))
            clusters[j].append(t)
            lab_list.append(j)
        labels = np.asarray(lab_list)

        # --- Track how many labels changed ---
        if track_changes and prev_labels is not None:
            changes = int(np.sum(labels != prev_labels))
            changes_per_iter.append(changes)

        # --- Early stop: labels unchanged ---
        if prev_labels is not None and np.array_equal(labels, prev_labels):
            converged = True
            break
        prev_labels = labels

        # --- Update step: Fréchet mean per cluster ---
        new_centers = []
        for members in clusters:
            if members:
                try:
                    m = frechet_mean_jpype(
                        members,
                        num_iter=mean_iters,
                        cauchy_len=5,
                        epsilon=epsilon,
                        return_java=False
                    )
                    new_centers.append(clean_newick(m))
                    print(f"Entered Try Block{it}")
                except Exception as e:
                    
                    # fallback: if mean computation fails, pick random member
                    new_centers.append(random.choice(members))
    
                    if hasattr(e, 'message'):
                        print(e.message)
                        #print(members)
                    else:
                        print(e)
                        #print(members)
                    
            else:
                # handle empty cluster by reseeding
                new_centers.append(random.choice(trees))

        # Secondary stop: centers unchanged (string equality)
        if all(new_centers[i] == centers[i] for i in range(k)):
            centers = new_centers
            converged = True
            break

        centers = new_centers

    # --- Compute final inertia (sum of squared distances to closest center) ---
    inertia = 0.0
    for t in trees:
        dmin = min(bhv_safe(c, t) for c in centers)
        inertia += dmin * dmin

    return inertia, it, converged, changes_per_iter

# ====== Elbow up to k=MAX_K (with optional change tracking) ======
Ks = list(range(1, min(MAX_K, n-1) + 1))
inertias = []
changes_by_k = {}  # k -> list of "# of trees that switched clusters" per iteration

print(f"\n=== Elbow up to k={MAX_K} | epsilon={EPSILON} | mean_iters={MEAN_ITERS} | seed={SEED} ===")
for k in Ks:
    print(f"k={k}...")
    inertia_k, it_used, conv, changes = kmeans_frechet(
        trees,
        k,
        iters=ITERS,
        epsilon=EPSILON,
        seed=SEED,
        mean_iters=MEAN_ITERS,
        track_changes=TRACK_CHANGES   # <--- collect changes during the elbow run
    )
    print(f"  converged={conv}  iters_used={it_used}  inertia={inertia_k:.20f}")
    inertias.append(inertia_k)
    if TRACK_CHANGES:
        changes_by_k[k] = changes  # may be empty if it stabilized immediately

# elbow heuristic (farthest point)
p1 = np.array([Ks[0], inertias[0]], float)
p2 = np.array([Ks[-1], inertias[-1]], float)
v  = p2 - p1
dists = [abs(np.cross(v, np.array([k, s]) - p1)) / (np.linalg.norm(v)+1e-12)
         for k, s in zip(Ks, inertias)]
elbow_k = Ks[int(np.argmax(dists))]

print("\nElbow results:")
print("k tested :", Ks)
print("inertias :", [round(x, 20) for x in inertias])
print("elbow_k  :", elbow_k)

# Inertia curve
plt.figure()
plt.plot(Ks, inertias, marker="o")
plt.xlabel("k"); plt.ylabel("Inertia")
plt.title(f"Elbow (k≈{elbow_k}) — {PATH.split('/')[-1]} [n={n}]")
plt.grid(True, linestyle="--", alpha=0.3); plt.xticks(Ks); plt.tight_layout(); plt.show()

# Reassignment curves (only if enabled)
if TRACK_CHANGES:
    # one plot with a line per k
    any_lines = False
    plt.figure()
    for k in Ks:
        ch = changes_by_k.get(k, [])
        if ch:
            it_axis = list(range(2, 2 + len(ch)))  # first comparison happens at iter 2
            plt.plot(it_axis, ch, marker="o", label=f"k={k}")
            any_lines = True
        else:
            print(f"[k={k}] No changes recorded (possibly converged right after first assignment).")
    if any_lines:
        plt.xlabel("Iteration")
        plt.ylabel("# trees that switched clusters vs. previous iter")
        plt.title(f"Reassignments per iteration (eps={EPSILON}, seed={SEED})")
        plt.grid(True, linestyle="--", alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.show()
