import jpype

# === Config ===
STURM_MEAN_JAR = "SturmMean_201102.jar"

def start_jvm():
    """Start the JVM with SturmMean on the classpath (no-op if already started)."""
    if not jpype.isJVMStarted():
        jpype.startJVM(classpath=[STURM_MEAN_JAR])

def shutdown_jvm():
    """Shut down the JVM (safe no-op if not running)."""
    if jpype.isJVMStarted():
        jpype.shutdownJVM()

"""
    Compute the Fréchet mean of rooted Newick trees via JPype (SturmMean.jar).

    Parameters
    ----------
    trees : list[str]
        Rooted Newick strings. A trailing ';' will be added if missing.
    num_iter : int
        Maximum number of iterations (default 10000).
    cauchy_len : int
        Length of Cauchy sequence used for convergence (default 10).
    epsilon : float
        Convergence threshold (default 0.05). Larger = faster but looser.
    outfile : str or None
        Optional path to write the mean; None disables file output.
    display_iter : int
        Display current approximation every N iterations (0 disables).
    display_start : int
        Start displaying after this iteration (only if display_iter > 0).
    return_java : bool
        If True return the Java PhyloTree; otherwise return a Newick string.

    Returns
    -------
    Java PhyloTree or str
        Java object if return_java=True, otherwise a Newick string (ends with ';').
    """

def frechet_mean(trees, num_iter=10000, cauchy_len=10, epsilon=0.05,
                       outfile=None, display_iter=0, display_start=1,
                       return_java=False):
    
    start_jvm()

    PhyloTree = jpype.JClass("distanceAlg1.PhyloTree")
    CentroidMain = jpype.JClass("centroid.CentroidMain")

    # Normalize semicolons and convert to Java array
    java_trees = jpype.JArray(PhyloTree)(
        [PhyloTree(t if t.endswith(";") else t + ";", True) for t in trees]
    )

    result = CentroidMain.getCentroidViaRandPermCauchy(
        java_trees, num_iter, cauchy_len, float(epsilon),
        outfile, display_iter, display_start
    )

    if return_java:
        return result

    s = str(result.toString()).strip()
    if not s.endswith(";"):
        s += ";"
    return s