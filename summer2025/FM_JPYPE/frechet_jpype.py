import jpype
from pathlib import Path

# === Config ===
STURM_MEAN_JAR = str(Path("SturmMean_201102.jar").resolve())
GTP_JAR        = str(Path("gtp.jar").resolve())


"""
    Start the Java Virtual Machine (JVM) with both JARs on the classpath.
    - This uses JPype to bridge Python and Java.
    - If the JVM is already running, this does nothing (no-op).
    """
def start_jvm():
    if not jpype.isJVMStarted():
        jpype.startJVM(classpath=[STURM_MEAN_JAR, GTP_JAR])
"""
    Shut down the JVM if it is currently running.
    - Useful for cleanup between runs or in interactive environments.
    """
def shutdown_jvm():
    if jpype.isJVMStarted():
        jpype.shutdownJVM()
"""
    Compute the Fréchet mean (centroid) of a set of phylogenetic trees
    using the Java implementation of randomized Cauchy sampling.
    
    Parameters
    ----------
    trees : list[str]
        A list of tree strings in Newick format.
    num_iter : int
        Number of optimization iterations (controls accuracy/runtime).
    cauchy_len : int
        Length of the Cauchy sequence used for centroid estimation.
    epsilon : float
        Convergence tolerance. Smaller values => more precise, slower.
    outfile : str or None
        Optional file path to save intermediate results (Java-side).
    display_iter : int
        How often (in iterations) to print progress. 0 disables output.
    display_start : int
        Iteration number at which display starts.
    return_java : bool
        If True, return the raw Java object instead of a Newick string.

    Returns
    -------
    str or JavaObject
        - Newick string of the Fréchet mean tree (default)
        - Raw Java object if `return_java=True`.
    """
def frechet_mean_jpype(trees, num_iter=10000, cauchy_len=10, epsilon=0.05,
                       outfile=None, display_iter=0, display_start=1,
                       return_java=False):
    start_jvm()

    PhyloTree     = jpype.JClass("distanceAlg1.PhyloTree")
    CentroidMain  = jpype.JClass("centroid.CentroidMain")
    

    java_trees = jpype.JArray(PhyloTree)(
        [PhyloTree(t if t.endswith(";") else t + ";", True) for t in trees]
    )

    result = CentroidMain.getCentroidViaRandPermCauchy(
        java_trees, int(num_iter), int(cauchy_len), float(epsilon),
        outfile, int(display_iter), int(display_start)
    )

    if return_java:
        return result

    s = str(result.getNewick(True))
   
    return s
