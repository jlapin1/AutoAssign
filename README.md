This project, completed between fall 2017 - spring 2018, seeked to automate the spectral assignment process for solid state protein NMR.
The result of the project was the following paper, published in JMR

Lapin, J., Nevzorov, A.A. (2018) Automated assignment of NMR spectra of macroscopically oriented proteins using simulated annealing. Journal of Magnetic Resonance, 293, 104-114

Following publication, the associated computer codes were turned into useable software, operable from a linux shell.
The code is as following:

1. speced.py
Python script that uploads NMR spectra and allows user to interactively edit the spectrum so that the automated search can work correctly.

2. mc.c
C source code to run Monte Carlo/Simulated annealing simulation as fast as possible on the edited spectra. Must be compiled into executeable to run.

3. anal.py
Python script that allows user to interactively analyze and filter results from the output of the C program
