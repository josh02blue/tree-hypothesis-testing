# Minimal example: compute a Fréchet mean from two rooted trees using JPype.
from frechet_jpype import frechet_mean, shutdown_jvm

trees = [
    "((A:0.1,B:0.2):0.3,(C:0.2,D:0.1):0.3);",
    "((A:0.12,B:0.18):0.28,(C:0.22,D:0.09):0.31);",
]

mean_newick = frechet_mean(
    trees,
    # defaults shown; tweak for speed/accuracy
    num_iter=10000,  
    cauchy_len=10,
    epsilon=0.05,
    outfile=None,
    display_iter=0,
    display_start=1,
    return_java=False
)

print("Fréchet mean (Newick):")
print(mean_newick)

# Optional cleanup (only when your whole program is done)
shutdown_jvm()