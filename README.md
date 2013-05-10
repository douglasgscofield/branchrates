branchrates
===========

Branchrates calculates branch-specific evolutionary rates given a set of taxa, a phylogenetic
tree for those taxa, and a set of homologous traits across all taxa.  It has been quiet for a while.

Load input data by calling functions, for taxa
(`taxon.matrix.read_file()`), a phylogenetic tree for
those taxa with branch lengths (`tree.read_internal_tree`), a set of homologous
traits across all given taxa (trait_matrix.read_file), and a mapping of
branch-level parameters to true parameters for each branch of the tree
(`ratep_map.read_file()`).

Please read the PDF documentation for an overview of mapping concepts, as my
use of that introduces probably the most confusing aspects of the
implementation.  Each branch of the tree has two parameters for character
evolution, a forward rate (character gain) and a backward rate (character
loss).  The pair of rates are managed as a unit, with each branch having an id
for the pair and a pointer to the single instantiation of `BranchRateManager`
which is queried for the rates.  `BranchRateManager` manages the mapping between
these branch-level parameters and the "true" parameters.  The "true" parameters
are kept in a private instantiation of `RatePVector`, which is initialized based
on the mapping described in a `RatePMap` instantiation, by calling
`BranchRateManager.allocate_ratep_from_map()`.

When computing the likelihood of a tree, the branch-level rates are retrieved
from `BranchRateManager` using the id for the pair.  In contrast, parameter
adjustments made while maximizing the likelihood are done directly to the
`RatePVector` through the `BranchRateManager`.  Either way, the value is ultimately
held in the `RatePVector`, it's just that the branch-level rates are abstracted
via the mapping.

The only current maximization implementation that I have confidence in for a
large number of parameters is the Nelder-Mead downhill simplex method, as
implemented in the `ML_multi_DownhillSimplex` class.  The current implementaton
is not licences for public distribution and will be replaced in favor of the one 
from the Gnu Scientific Library.


Current Setup
-------------

The current branchrates uses a dataset from the Koonin lab to calculate a
number of evolutionary parameters for intron birth and death.  The phylogenetic
tree used is one I determined using PAUP, with branch lengths multiplied by
100.  This tree is kept in an internal array because I don't yet have a method
for reading an external tree, one is partially present but not completed.  The
list of taxa, the trait matrix, and the parameter mapping are all kept in
external files as is clear from `main()`.  Output generated includes a summary
of the `TraitMatrix`, the `RatePMap`, the `PhyloTree`, and then the output of
the maximization process.  Each "`amoeba()`" line includes the low and high
likelihoods computed among the vertices of the simplex, which can be used to
follow the minimization process.  I've allowed for 30 restarts, each restart is
initiated when the difference between the likelihoods within the simplex falls
below a minimum.  At the end of each restart, the current estimates of the
parameter values are printed.  This takes quite a while to finish, given the 30
restarts.  The method `ML_multi::profile_likelihoods_print` prints out
parameter profiles surrounding each parameter estimate based on chi-square
likelihood ratio criteria.

