import jpype
import jpype.imports
import time
import random
import statistics
import subprocess
from pathlib import Path

# === Config ===
GTP_JAR = "gtp.jar"                 # path to gtp.jar
INPUT_FILE = "rootedtree_1.txt"     #file with ~1000 rooted Newick trees, one per line
LIMIT_TREES = 1000                  #read at most this many lines/trees from the file
PAIRS_MODE = "random"             #"adjacent" or "random"
N_PAIRS = 500                       #how many pairs to time (<= number of available pairs)
ROOTED = True                       #your trees are rooted
DIST_TOL = 1e-9                     #equality tolerance for distance agreement check
SEED = 42                           #for reproducible random pairs

# --- JPype method ---
def bhv_jpype(t1: str, t2: str, rooted: bool = True) -> float:
    
    if not jpype.isJVMStarted():
        jpype.startJVM(classpath=[GTP_JAR])

    PolyMain = jpype.JClass("polyAlg.PolyMain")
    PhyloTree = jpype.JClass("distanceAlg1.PhyloTree")

    t1 = t1.strip()
    t2 = t2.strip()
    if not t1.endswith(";"): t1 += ";"
    if not t2.endswith(";"): t2 += ";"

    pt1 = PhyloTree(t1, rooted)
    pt2 = PhyloTree(t2, rooted)
    geo = PolyMain.getGeodesic(pt1, pt2, None)
    return float(geo.getDist())

# --- CLI method ---
def bhv_cli(t1: str, t2: str) -> float:
    t1 = t1.strip()
    t2 = t2.strip()
    if not t1.endswith(";"): t1 += ";"
    if not t2.endswith(";"): t2 += ";"

    # Pass Newick strings directly (no temp files)
    cmd = ["java", "-jar", GTP_JAR, "-t", t1, "-o", "out.txt", t2]
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return float(res.stdout.strip())

def load_trees(path: str, limit: int):
    """Read up to `limit` trees (non-empty lines)."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    trees = []
    with p.open() as f:
        for line in f:
            ln = line.strip()
            if ln:
                trees.append(ln)
                if len(trees) >= limit:
                    break
    if len(trees) < 2:
        raise ValueError("Need at least 2 trees to make a pair.")
    return trees

def make_pairs(trees, mode: str, n_pairs: int, seed: int = 42):
    """Return a list of (t1, t2) pairs."""
    pairs = []
    if mode.lower() == "adjacent":
        # consecutive pairs: (0,1), (1,2), ...
        for i in range(len(trees) - 1):
            pairs.append((trees[i], trees[i + 1]))
    elif mode.lower() == "random":
        random.seed(seed)
        idxs = list(range(len(trees)))
        for _ in range(len(trees) * 2):  # oversample pool
            i, j = random.sample(idxs, 2)
            pairs.append((trees[i], trees[j]))
    else:
        raise ValueError("PAIRS_MODE must be 'adjacent' or 'random'.")

    return pairs[:min(n_pairs, len(pairs))]

def fmt_stats(arr):
    if not arr: return "n/a"
    return f"mean={statistics.mean(arr):.6f}s, median={statistics.median(arr):.6f}s, min={min(arr):.6f}s, max={max(arr):.6f}s"

def main():
    # Load trees
    trees = load_trees(INPUT_FILE, LIMIT_TREES)
    print(f"Loaded {len(trees)} trees from {INPUT_FILE}")

    # Build pairs
    pairs = make_pairs(trees, PAIRS_MODE, N_PAIRS, SEED)
    print(f"Timing {len(pairs)} pairs (mode={PAIRS_MODE})")

    # Warm-up JPype once (don’t count JVM startup in timings)
    jpype_ready = False
    try:
        _ = bhv_jpype(pairs[0][0], pairs[0][1], ROOTED)
        jpype_ready = True
    except Exception as e:
        print(f"JPype warm-up failed: {e}")

    # Time JPype
    jp_dists, jp_times = [], []
    if jpype_ready:
        for a, b in pairs:
            t0 = time.perf_counter()
            d = bhv_jpype(a, b, ROOTED)
            t1 = time.perf_counter()
            jp_dists.append(d)
            jp_times.append(t1 - t0)

    # Time CLI
    cli_dists, cli_times = [], []
    for a, b in pairs:
        t0 = time.perf_counter()
        d = bhv_cli(a, b)
        t1 = time.perf_counter()
        cli_dists.append(d)
        cli_times.append(t1 - t0)

    # Compare distances
    if jpype_ready:
        agree = all(abs(x - y) < DIST_TOL for x, y in zip(jp_dists, cli_dists))
        print(f"\nDistance equality (JPype vs CLI): {agree}")
        if not agree:
            for i, (x, y) in enumerate(zip(jp_dists, cli_dists), 1):
                if abs(x - y) >= DIST_TOL:
                    print(f"  mismatch #{i}: JPype={x}  CLI={y}")

    # Timing summary
    print("\n=== Runtime Summary ===")
    if jpype_ready:
        print(f"JPype  : {fmt_stats(jp_times)}  (n={len(jp_times)})")
    else:
        print("JPype  : failed to initialize; no timings.")
    print(f"CLI    : {fmt_stats(cli_times)}  (n={len(cli_times)})")

    # Clean JVM
    if jpype_ready:
        import jpype
        if jpype.isJVMStarted():
            jpype.shutdownJVM()

if __name__ == "__main__":
    main()