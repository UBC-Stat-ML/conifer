TODO
----

### Question: to normalize or not to normalize the rate matrix

**OK:** Are there advantages of a non-identifiable parameterization? E.g. is impute the path, could do analytic posterior with a dirichlet-gamma parameterization. **NO:** might be overly complicated. Better stick with the current way.

### Clean implementation of CTMC

**DONE:** Start with CTMC.java. Use EJML and diagonalization. Very basic rate matrix (fixed) for now (already done in scaffold?). 

### Test phylogenetics code

1. Write forward simulation code (done for prior, just need it for likelihood) (**DONE**). 
2. Make TreeLikelihood implement GenerativeFactor. (**DONE**)
3. Change name of TreeLikelihood to UnrootedTreeLikelihood (**DONE**)
4. Write phylo test case
  - Some basic processors on trees (**DONE**)
  - That's it! Run the test! (**DONE**)
  - Check what went wrong. (**DONE**)
5. A few more stats  <-- entropy on clade distribution? pr of largest non-trivial clade?
**OK, one more implemented, but how to ensure our test is non-trivial?**

### Investigate gradient code

See why penalized likelihood is decreasing

### Compare to MrBayes by using logGamma priors?

### Import HMC project

Test on simple example.

### Integerate gradient code

Start with naive (slow) HMC way

### Write fast uniformization-based sampler

### Build consensus code

Existing code bases:

- Beast: seems to require rooting
- Mesquite: seems complex to interface
- Hashcs: needs to be installed, does not report clade support
- Sumt: needs to be installed, depends on having a nexus format
- Legacy: quadratic time, tied to old datastructures

Write new: this is done in two stages:

1. collect clades statistics
2. build consensus

For 1, can use the sum product architecture (as was done for fast-tree-metric in legacy). In addition, assign random long to each leave, and have the hash of internal be the sum of the two longs. Put that in hash. Do it twice with different random longs at leaves. Pr of no collision is greater than (1-np)^n where p is 2^-64 and n is the number of nodes in the tree. So this should work with very high pr even for large tree. 

For 2, pick arbitrary leaf v and orient each tree from it. Traverse the tree and record parent child relationships of consensus nodes only in one shared network. Assign max distance to v to each node in the consensus network. Use this to remove redundant edges. This should yield a tree (prove it).

### Fix real variable sampler

- Back to std one
- Do something better?

### Add list of Factor syntax

### Better default processor

E.g. variables should produce graphs

### Check points, initializations

### Command line interface

### Remove proposal


### Make creation better

- initialization in test cases a bit ugly

### Fix tutorial

Try to get students involved!

### Add other distributions

Try to get students involved!

