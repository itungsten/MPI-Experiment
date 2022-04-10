# Parallel-and-Distributed-Computing-MPI-Experiment

##  Introduction
*Parallel and Distributed Computing*  of UESTC

Course Experiment: Implements of Sieve of Eratosthenes, with optimizer and improvement

## Usage

1. make
2. mpiexec -n \<num_threads\> ./optimizer1 \<limit\>  
3. mpiexec -n \<num_threads\> ./optimizer2 \<limit\> 
4. mpiexec -n \<num_threads\> ./optimizer3 \<limit\>  [chunk_size]
5. mpiexec -n \<num_threads\> ./optimizer4 \<limit\>  [chunk_size]
6. make clean

e.g.

```bash
make
mpiexec -n 4 ./optimizer1 1000000000
mpiexec -n 4 ./optimizer4 1000000000
mpiexec -n 4 ./optimizer4 1000000000 23333
make clean
```

## Environment

Base on 
* Microsoft Windows 10 Pro 10.0.19044 Build 19044
* WSL Ubuntu 18.04 bionic
* Intel Core i7-8565U @ 8x 1.992GHz (4 physical cores)
* MPICH 3.3a2

## File Description
* base.cpp: The origin source code provided(corrected some bugs),
* optimizer1.cpp: base + 2's reduced,
* optimizer2.cpp: opt1 + broadcast reduced,
* optimizer3.cpp: opt2 + cache reduced,
* optimizer4.cpp: opt3 + 3's and 5's reduced + union + O3.
