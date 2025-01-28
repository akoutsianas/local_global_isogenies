
## Local-global principle for isogenies of elliptic curves over quadratic fields

  

This repository accompanies the paper *`Local-global principle for isogenies of elliptic curves over quadratic fields'* by Stevan Gajovic, Jeroen Hanselman and Angelos Koutsianas. (arXiv link will follow).


### Contents of the repository

We give a short overview of the files in our repository:


- **XD10_mordell_weil_torsion.m** (Section 4 of the paper) In this file we determine a finite index subgroup of the Mordell-Weil group of the Jacobian of $X_{D_{10}} (11)$.

- **models_and_maps_computations.m** (Section 3 of the paper)  We compute a model for $X_{D_{10}}(11)$ using Assaf's code and show that it is isomorphic to a simpler model. We furthermore compute the Atkin-Lehner involutions $w_{11}$, $w_{121}$ and the quotient maps to $X_0^+(121)$ and $X^+_{D_{10}}(11)$. 

- **mordell_weil_quotient.m** (Section 4 of the paper) In this file we compute the Mordell-Weil group and the rational points on $X^+_{D_{10}}(11)$.

- **symmetric_chabauty.m** This file contains an implementation of the relative symmetric Chabauty method and the Mordell-Weil sieve for the curve $X_{D_{10}}(11)$.

- **symmetric_chabauty_computations.m** (Section 7 of the paper) In this file we combine everything and apply the relative symmetric Chabauty method  implemented in symmetric_chabauty.m and the Mordell-Weil sieve to $X_{D_{10}}(11)$.

All of our computations were done on a machine running Ubuntu 22.04.1 with an Intel i7-7700 processor, 4 cores @ 3.60GHz, 32GB RAM.

## Installation instructions

First make sure you have the latest version of Magma installed. (The code was written for Magma V2.28-16).

Then use the command

```
git clone https://github.com/akoutsianas/local_global_isogenies.git
```
in the folder you want to download the files to.

### Section 3 of the paper
For the computations we did in this section we used code written by Eran Assaf. In order to install the ModFrmGL2 package by Eran Assaf use the command:
```
git clone https://github.com/assaferan/ModFrmGL2
```
In the following we will assume that this was done in the same folder where the local_global_isogenies folder is located. After this perform the commands:
```
git submodule init
git submodule update
```
This will finish the installation of Assaf's package.
Now go to the **local_global_isogenies** folder and start Magma. After this we attach the spec of Assaf's code. This can be done by using the command:
```
AttachSpec("../ModFrmGL2/ModFrmGL2/ModFrmGL2.spec");
```
(P.S. Depending on your setup the path to the spec might need to be modified.)
After this you can run the file by either using:
```
iload "models_and_maps_computations.m";
```
or copying the contents of the file line by line if you are only interested in parts of the computation.

### Section 4 of the paper
Computing the Torsion part is done by going to the **local_global_isogenies** folder, starting Magma and by calling for example:
```
iload "XD10_mordell_weil_torsion.m";
```
or copying the contents of the file line by line if you are only interested in parts of the computation.

For the rest of the computations in this section we need some files written by Samir Siksek. Download **add.m** and **g2-jac.m** from https://github.com/samirsiksek/siksek.github.io/tree/main/progs/chabnf and place them in the **local_global_isogenies** folder.

Now go to the **local_global_isogenies** folder and start Magma. After this we attach Siksek's code. This can be done by using the command:
```
Attach("add.m");
Attach("g2-jac.m");
```
(P.S. Depending on your setup the path to the spec might need to be modified.)

After this you can run the file by either using:
```
iload "mordell_weil_quotient.m";
```
or copying the contents of the file line by line if you are only interested in parts of the computation.

The computations took a few seconds on our system.

### Section 7 of the paper

Performing the symmetric Chabauty computation is done by going to the **local_global_isogenies** folder, starting Magma and calling:
```
Attach("symmetric_chabauty.m");
iload "symmetric_chabauty_computations.m";
```
or copying the contents of the file line by line if you are only interested in parts of the computation.

(P.S. Depending on your setup the path to the spec might need to be modified.)
The computations took about two minutes on our system.
