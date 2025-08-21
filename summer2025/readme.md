- SturmMean : java jar file that contains code for computing (Frechet) means of a set of trees
<<<<<<< HEAD
- trees_to_edge_lengths.py  : Python script doing multiple things, including simulating gene trees for a given species tree and computing the Frechet mean of the gene trees

How to run Permutation Test

- Generate Frechet Means and Gene Trees
In trees_to_edge_lengths.py, generate:
The gene tree files for two different species
The corresponding Frechet mean files for each species

- Update Paths in permutation_test_gene_trees.py
Edit the user settings section at the top of permutation_test_gene_trees.py to include the correct relative paths for:

- The gene tree files (SPECIES1_TREES, SPECIES2_TREES)

- The precomputed Frechet mean files (REAL_MEAN_1, REAL_MEAN_2)

- Run the Permutation Test File
=======
- trees_to_edge_lengths.py : Python script doing multiple things, including simulating gene trees for a given species tree and computing the Frechet mean of the gene trees

-How to run Permutation Test

Generate Frechet Means and Gene Trees
In trees_to_edge_lengths.py, generate:

The gene tree files for two different species

The corresponding Frechet mean files for each species

Update Paths in permutation_test_gene_trees.py
Edit the user settings section at the top of permutation_test_gene_trees.py to include the correct relative paths for:

The gene tree files (SPECIES1_TREES, SPECIES2_TREES)

The precomputed Frechet mean files (REAL_MEAN_1, REAL_MEAN_2)

Run the Permutation Test File
>>>>>>> beb4922 (add bhv folder)
