# A MATLAB implementation of greedy algorithms for unsupervised variable selection
The repository contains the source code of the paper: Federico Zocco, Marco Maggipinto, Gian Antonio Susto and Seán McLoone, 2023. "Lazy FSCA for unsupervised variable selection", Engineering Applications of Artificial Intelligence, vol. 124, 106624. 

If you use this code, please consider citing the corresponding paper mentioned above.

# Dependencies of folders and files
First, the user should refer to Simulations.m as it demonstrates the use of the code and generates the paper results. The seven greedy algorithms covered in the source paper are: (1) fsca.m, (2) lazy_fsca.m, (3) OPFS.m, (4) fosmod.m, (5) GPR.m, (6) fsca_fsfp.m and (7) ufs.m. 

The datasets "Sim1" and "Sim2" are generated by "DataGenerator.m". "PitProps.mat" is uploaded on this repository. The datasets "Semiconductor" and "MofERdata" (called "PlasmaEtch" in the source paper) are not public for confidentiality reasons. The remaining datasets are not uploaded on this repository due to memory constraints. The user can download them from the links provided in the source paper. 

