# Introduction

This repository extends the Hatree Fock method programme based on the programme provided in the Modern Quantum Chemistry by Attila Szabo.
The origin code in Szabo's book performs single zeta -NG(N ranges from 1 to 3) for H-He^{+}, which contains 2 electrons.

In this repository, specifically main.f90, I extended the code in Szabo's book, making it capable of performing arbitrary number of
atoms and electrons. 
The user only have to define certrain constant variables and provide a formatted file indicating the configuration of atoms.

# Usage of programme 

As the user can see in main.f90, the following variables should be provided

num\_basis: total number of basis (if num\_basis .neq. num\_atom,
the user must change the init()function to generate other basis functions)

num\_occupied\_orb: number of occupied orbitals / 2 because of closed-shell

num\_atom : total number of atoms

max\_iter : maximum of iterations

Please don't forget to change the path of file containing the configuration.
The format of the file can be found in geomtry\_rep.xyz, which is the standard
import form of fhi-aims(an outstanding ab initio calculation software). 

After determining the necessary information mentioned above, the user just have to compile
the programme. For example:

```shell
ifort main.f90 -llapack -lblas -o myprogramme
```

And run the programme:

```shell
./myprogramme
```

Note that the results of the programme I uploaded is the same with that in Szabo's book, which verifies that 
there is highly likely no error in my programme.
# Notes

The default basis in the main.f90 is single zeta-3G 1s for H and He.
And the code deals with obviously the closed-shell problem.
If other elements is needed, the user should provide coeffcients and exponent by him/herself.
See subroutine init()

Please note that there is no garantee of convergence for complex systems, because the completeness of the
basis is inadiquate.
