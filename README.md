-------------------------------------------------------
Optimized Coordinate Descent Massive MU-MIMO simulator
-------------------------------------------------------
(c) 2016 Christoph Studer and Michael Wu 
e-mail: studer@cornell.edu 
-------------------------------------------------------

# Important information:

If you are using the simulator (or parts of it) for a publication, then you MUST cite our paper:

@article{WuDickCavallaroStuder:2016,
  author={M. Wu and C. Dick and J. R. Cavallaro and C. Studer},
  journal={IEEE Transactions on Circuits and Systems I},
  title={High-Throughput Data Detection for Massive MU-MIMO-OFDM using Coordinate Descent},
  year={2016}
}

and clearly mention this in your paper (also, you may want to add the month, page numbers, issue number, etc. in the above BibTex entry). More information about our research can be found at 

http://vip.ece.cornell.edu

# How to start a simulation:

Simply run 

>> OCD_sim

which starts a simulation in a 128 BS antenna, 8 user massive MU-MIMO system with 16-QAM using various detectors, including ZF and MMSE detection, the SIMO lower bound, and optimized coordinate descent with BOX and MMSE regularization (OCD performs only 2 iterations; note that the performance is almost the same as for the MMSE detector and the SIMO lower bound!). The simulator runs with predefined simulation parameters. You can provide your own system and simulation parameters by passing your own ?par?-structure (see the simulator for an example). 

We highly recommend you to execute the code step-by-step (using MATLAB?s debug mode) in order to get a detailed understanding of the simulator. 

# Version 0.1 (September 5, 2016) - studer@cornell - initial version for public access