# Minimal BHV example using Newick strings

from bhv_jpype import ensure_jvm, bhv_distance, shutdown_jvm

if __name__ == "__main__":
    ensure_jvm()  # starts JVM with gtp.jar

    t1 = "(A,(B,C));"     # rooted Newick
    t2 = "((A,B),C);"     # rooted Newick

    d = bhv_distance(t1, t2, rooted=True)
    print("BHV distance:", d)

    shutdown_jvm()        # optional