#DNA-Metropolis
This package provides a Monte Carlo simulation of the 3D DNA formation confined within the nucleus based on the random loop model for long polymers by Bohn et al. [1], which has already been successfully used to describe folding of chromatin in human cells [2].

1. Bohn, M., Heermann, D. W. and van Driel, R. (2007). Random loop model for long polymers. Phys Rev E Stat Nonlin Soft Matter Phys 76, 051805.

2. Mateos-Langerak, J., Bohn, M., de Leeuw, W., Giromus, O., Manders, E. M. M., Verschure, P. J., Indemans, M. H. G., Gierman, H. J., Heermann, D. W., van Driel, R. and Goetze, S. (2009). Spatially confined folding of chromatin in the interphase nucleus. Proc Natl Acad Sci U S A 106, 3812–3817.



##Compile and run
Make sure that g++ and boost are installed, on debian-like systems: 

```apt-get install build-essential libboost-all-dev)```

Then compile with

```g++ -std=c++0x -O3 -Wall main.cpp dna_metropolis.cpp -o dna_metropolis```

Finally, run the simulation with `./dna_metropolis`.

