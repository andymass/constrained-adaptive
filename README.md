# Constrained adaptive sensing code

Code for the pre-print [Constrained adaptive sensing](http://arxiv.org/abs/1506.05889) 
by Mark A. Davenport, Andrew K. Massimino, Deanna Needell, and Tina Woolf

BibTeX reference:
```
@article{DBLP:journals/corr/DavenportMNW15,
   author    = {Mark A. Davenport and
                Andrew K. Massimino and
                Deanna Needell and
                Tina Woolf},
  title     = {Constrained adaptive sensing},
  journal   = {CoRR},
  volume    = {abs/1506.05889},
  year      = {2015},
  url       = {http://arxiv.org/abs/1506.05889},
  timestamp = {Wed, 01 Jul 2015 15:10:24 +0200},
  biburl    = {http://dblp.uni-trier.de/rec/bib/journals/corr/DavenportMNW15},
  bibsource = {dblp computer science bibliography, http://dblp.org}
}
```

## Requirements

- SPGL1
    https://www.math.ucdavis.edu/~mpf/spgl1/

Van Den Berg, E., & Friedlander, M. P. (2007). SPGL1: A solver for large-scale
sparse reconstruction.

Van Den Berg, E., & Friedlander, M. P. (2008). Probing the Pareto frontier for
basis pursuit solutions. SIAM Journal on Scientific Computing, 31(2), 890-912.

- TFOCS
    http://cvxr.com/tfocs/

Becker, S. R., Candès, E. J., & Grant, M. C. (2011). Templates for convex cone
problems with applications to sparse signal recovery. Mathematical programming
computation, 3(3), 165-218.

- Rice Wavelet toolbox (RWT) 
    http://dsp.rice.edu/software/rice-wavelet-toolbox

- Image Processing Toolbox

- `brain.mat' data file for two-dimensional experiments
    http://www.eecs.berkeley.edu/~mlustig/CS/brain.mat
    http://www.eecs.berkeley.edu/~mlustig/CS.html

## Experiments

### 1 dimensional experiments

See `tester_UnifAndVDS.m` and `tester_n_UnifAndVDS.m`.

### 2 dimensional experiments
    
See `brain_tester.m'.

