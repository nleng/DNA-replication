#Replication
This program simulates the time evolution of DNA-replication on a 1D DNA-string. To generate a realistic 3D-Simulation is has to be combined with the DNA-Metropolis package.

##Compile and run

Make sure that the [simtools library](../simtools/) is properly installed.

Compile:

```g++ -O3 -std=c++11 annihilation.cpp boundary.cpp exceptions.cpp fork.cpp parameter_set.cpp replicator.cpp twoe.cpp main.cpp -o DNA_replication -lsimtools```

Run with: 
```
./DNA_replication
```



