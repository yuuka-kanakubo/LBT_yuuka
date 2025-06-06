# LBT (Linear Boltzmann Transport) WARNING: Unofficial version!! Please let me know if you want to run this code.

A LBT Monte Carlo model describes jet propagation and interaction with the quark-gluom plasma (QGP) in relativistic heavy-ion collisions.
For the latest full paper about LBT, see Tan Luo(CCNU, Wuhan, Inst. Part. Phys. and Santiago de Compostela U., IGFAE), Yayun He(CCNU, Wuhan, Inst. Part. Phys. and South China Normal U.), Shanshan Cao(Shandong U.), Xin-Nian Wang(CCNU, Wuhan, Inst. Part. Phys. and LBL, Berkeley), Phys.Rev.C 109 (2024) 3, 034919.

This version is refactored version of the LBT by Yuuka Kanakubo (RIKEN iTHEMS, LBL, UC Berkeley).  


- Build with CMake
```
mkdir build
cd build
cmake ..
make
mv src/LBT ../
cd ..
```


- Run unit test (more to be added)
```
cd src/test
./RunTest.sh
```

