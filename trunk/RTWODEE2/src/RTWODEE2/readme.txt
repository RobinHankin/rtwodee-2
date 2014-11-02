If you already have compiled code then R-TWODEE2 simply requires Java 1.7 or above to run.
However to first compile the code you will need to have the JDK installed which can be downloaded from http://www.oracle.com/technetwork/java/javase/downloads/index.html

To check if you have JDK installed, type javac -version into the command line (must be version 1.7 or above).

Heavy gas dispersion using the shallow water approximations, project based on earlier TWODEE and TWODEE-2 software.

To run in command line:
1) cd to the directory the RTWODEE2 folder is in (NB: if in the Netbeans Project folders, this is located in the src directory)
2) Compile all classes with: javac RTWODEE2/*.java
3) Run code with: java RTWODEE2/TWODEE
4) Re-compiling with step 2 is only needed if the source code is changed, not if any of the data input files are changed

Editing a Simulation:
- The main settings of the TWODEE simulation is stored in dataFiles/genSettings.in
- The output for the simulation will appear in scenario folder selected once the simulation has been run (NB: if running in Netbeans, this will be build folder not the src folder)

Editing the Code:
- The Simulation class contains the calculations the simulation uses
- The runSimulation method contains the method calls to each calculation, these can be commented out or re-ordered if need be

Contacting Dr. Robin Hankin:
Email: hankin.robin@gmail.com

This project is available on Google Code at https://code.google.com/p/rtwodee-2/


