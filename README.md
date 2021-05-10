# rheoInterEHDFoam
Extention for rheoInterFoam to solve electric equations for two-phase visco-elastic flows

This solver has been developed as a part of my undergraduate thesis at Birla Institute of Technology and Science (BITS), Pilani. The work would not have been complete without the help of my supervisor Dr Shyam Sunder Yadav, Assistant Professor, Mechanical Engineering Department, BITS Pilani; who guided me at every step of the thesis.

**Prerequisites:**
Before installing the solver you need to install OpenFoam v7 from https://openfoam.org/version/7/ and rheoTool for OpenFoam v7 from https://github.com/fppimenta/rheoTool. 

**How to install rheoInterEHDFoam:**
1. Paste the rheoInterEHDFoam folder in rheoTool/of70/src/solvers folder.
2. Open the folder rheoInterEHDFoam in terminal.
3. Source OpenFoam and run _wmake_ command.
4. The solver is installed and ready to used.

The solver comes with a test case dropCoalescing.

**How to run test case:**
Go to the test case directory and open the terminal there. Run the following commands sequentially in the terminal:

1. blockMesh
2. setFields
3. rheoInterEHDFoam

The results can be viewed in ParaView, a software package shipped with OpenFoam to visulise the data.
