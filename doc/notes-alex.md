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
**OK for now, NNI code is pretty simple**

### Investigate gradient code

**Problem potentially identified:** normalization incorrect in the multi-category case? 

### Planning

Potential solution for above 

1. Fix normalization
   - Java **Too (programming) time consuming?**
   - Stan **Too much plumbing**
   - Just rescale the holding time statistics **Yes!**
2. Use only unnormalized matrices **No, too non-standard**

So to summarize, in the CTMC learning code, we will still assume a single rate matrix, but with two tricks:

- an actual category will be sampled for each site, so only the counts for this category will be accumulated
- waiting times can then be interprete rescaled with respect to the selected category (so that the full problem is seen as a single normalized matrix from the point of view of the CTMC training code)

**Actually:** wait for talking with Crystal, but there is probably no problem with the code on our side. 

**Cancel wait:** problem was probably that the rate over categories was not properly normalized there, making the true parameter not accessible in simulations.

### Implementing and testing new tree likelihood class

Next: 
- fixed rate matrix:
  - load and save
  - use easy kimura
  - std gamma rate cat **DONE**

### Simplest weight-based sampler

Use MH **DONE** simple tests needed on the basic infrastructure pre-HMC.

### Compare to MrBayes by using logGamma priors?

**NO, USED EXISTING TESTS INSTEAD**

### Import HMC project

Test on simple example. **DONE**

### Integerate gradient code

Start with naive (slow) HMC way **SKIPPED**

### Write fast uniformization-based sampler

**DONE**

---

  
  OK this test works now it seems.
  Next step:
     - before doing anything big: a bit of doc, version control, push maven
     - improvements in Blang?
     - some processor stuff to see if sampled matrices are reasonable in real data (and how fast? bunch of slow things in the expfam model that perhaps should be cached)
     - put back the modular part where json files can be used to specify features
     - document
     - rush final set of xp for the HMC stuff (cf MH, etc)
     - put back fancier tree moves already developed
     - need more processors (viz paths, params, ..)
     - par. tempering
     
---

### Refactor FactorGraph

**NO:** actually, the DiscreteFactorGraph does have some fields, i.e. keeping track of the number of sites. Also make add binary factor calls simplers.
 
- should be class instead of interface
- add invert in binary factors instead
- basically, main behavior is to pointwise multiply things
- remove newFactorGraph in interface EvolutionaryModel. called by sumprod


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


Topics of each paper
--------------------

### Models and algorithms for large substitution matrix estimation

**Journal:** Sys bio

**Topics:**

- Normalized exponential family model
- Applications: 
  - hierarchical, site-specific iid model
  - non-hierarchical: codon interaction evolution: features:
    - protein marginal features
    - codon marginal features
    - same for interaction features
  - non-exponential waiting times?
- Estimation
  - Frequentist: 
    - EM
    - gradient method
  - Bayesian:
    - HMC background; including new method based on Bayesian optimization
    - Auxiliary variable construction (applies Neal's paper)
    - Trick
    
    
Current bugs:

- behavior not deterministic

expected counts not the same across execution (but same within execution)
  suggests the input of the expected counts is diff; could be:
    tree
      different initially, resampled non-deterministically in between? 
        that seems to be it!
        OK: first problem identified: order of samplers
          Bug fixed: tracked down all newHashSet() and replaced by newLinkedHashSet()
    other?

- numerical error in ExpFamParameters
  fixed a NaN*0 = NaN -> NaN*0 = 0