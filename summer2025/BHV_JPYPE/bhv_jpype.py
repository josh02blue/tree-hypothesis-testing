import jpype, jpype.imports

# path to jar file
GTP_JAR = "gtp.jar"  


"""
Start the JVM with gtp.jar on the classpath if it is not already running.

Notes
-----
- The JVM can be started only once per Python process. Subsequent calls
    are safe: they will be no-ops.
- If you also need other JARs (e.g., SturmMean.jar), put *all* of them
    in the very first `jpype.startJVM(classpath=[...])` call in your script.
"""
def ensure_jvm():
    if not jpype.isJVMStarted():
        jpype.startJVM(classpath=[GTP_JAR])
"""
Shut down the JVM if running.

Notes
-----
Optional. You can omit this at program end.

"""
def shutdown_jvm():
    if jpype.isJVMStarted():
        jpype.shutdownJVM()
"""
Compute the BHV geodesic distance between two trees using gtp.jar.

Parameters
----------
tree1, tree2 : distanceAlg1.PhyloTree or str
    Either Java `PhyloTree` objects (recommended if you already have them),
    or rooted Newick strings. If strings are provided and the trailing ';'
    is missing, it will be added automatically.
rooted : bool
    Whether to parse strings as rooted trees (only used for string inputs).

Returns
-------
float
    The BHV geodesic distance.

Raises
------
Java exceptions (via JPype) if input Newick strings are malformed.
"""
def bhv_distance(tree1, tree2, rooted=True):
    
    ensure_jvm()
    PolyMain = jpype.JClass("polyAlg.PolyMain")
    PhyloTree = jpype.JClass("distanceAlg1.PhyloTree")

    # If both inputs are Java objects, skip parsing
    if hasattr(tree1, "getClass") and hasattr(tree2, "getClass"):
        pt1, pt2 = tree1, tree2
    else:
        t1 = (tree1 or "").strip()
        t2 = (tree2 or "").strip()
        if not t1.endswith(";"): t1 += ";"
        if not t2.endswith(";"): t2 += ";"
        pt1 = PhyloTree(t1, rooted)
        pt2 = PhyloTree(t2, rooted)

    geo = PolyMain.getGeodesic(pt1, pt2, None)
    return float(geo.getDist())

