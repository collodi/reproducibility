
Source files:
	
	HH_Blanco2016v1_*.cpp (main program)
	 - 100 excitatory HH neurons, all to all connected

	iappDist.cpp -> Implementation of functions to setup the Iapp randon values for 100, 200 and 400 neurons. 
		      Also a function to write the Iapp values in a .txt file

	SpikeTrain.cpp -> Class implementation. Class in charge to manipulalate the Spike train of the network. 


include files:
	iappDist.h -> Declaration of functions to setup the Iapp randon values for 100, 200 and 400 neurons. 
		      Also a function to write the Iapp values in a .txt file

	SpikeTrain.h -> Class declaration, which is a class in charge to manipulalate the Spike train of the network

 

Notes: 
 - You have to create the makefile
 - The code is ready to run for 8s of real simulation, dt = 0.01ms
 - In the folder in which the .exe is created must create a folder named "results"
 - Example how Call the program on the command line: 
	HH_Blanco2016v1_doubleP -case 21 -pEexcN 1 -nBurst 50 
 - In order to change the order of the values of Iapp, check the file iappDist.cpp


Software: 
 - Compiler: gcc version 6.3.0 (x86_64-posix-seh-rev1, Built by MinGW-W64 project)
 - boost_1_60_0
  
