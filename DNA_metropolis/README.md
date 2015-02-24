##Compile command

g++ -std=c++0x --shared -O3 -Wall -I/usr/include/python2.7 dna_metropolis.cpp get_plotpoints.cpp -lboost_python -fPIC -lpython2.7 -o libreplihelpers.so
