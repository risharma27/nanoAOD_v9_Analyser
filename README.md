## nanoAOD_v9_Analyser
### Basic analyser codes for nanoAOD root files compatible with nanoAOD version 9.

- The input ntuples to be analysed are stored in the ``nanoAOD_v9_Analyser/inputs`` directory.
- The analysis codes are stored in the ``nanoAOD_v9_Analyser/`` directory. The main code (which does analysis and fills histograms) is in the nano9Ana class. This class is declared and described in ``nano9Ana.h`` and ``nano9Ana.C``.

#### How to run the scripts

Open root prompt in this directory and load the script as follows:

```
[] .L nano9Ana.C+
```

This compiles the C code and makes it executable. Once loading is complete, execute the following in the ROOT prompt:

```
[] .x ana.C(1)
```
- ana.C is the driver script.
- The argument inside ana.C() stands for the different input files.
- The location to the input files should be given in the ana.C file before this step.
- This step generates a ``hst_<process name>.root`` file in ``nanoAOD_v9_Analyzer/hst_outputs`` directory and a ``sum_<process name>.txt`` file in ``nanoAOD_v9_Analyzer/sum_output`` directory.
