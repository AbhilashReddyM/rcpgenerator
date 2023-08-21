CPP = g++
FLAGS = -O3 -march=native

rcpgenerator.exe: rcpgenerator.o distance.o minimization.o histogram.o readwritedisplay.o twister.o
	$(CPP) $(FLAGS) rcpgenerator.o distance.o minimization.o histogram.o readwritedisplay.o twister.o -o rcpgenerator.exe

rcpgenerator.o: rcpgenerator.cpp
	$(CPP) $(FLAGS) -c rcpgenerator.cpp

twister.o: twister.cpp
	$(CPP) $(FLAGS) -c twister.cpp

histogram.o: histogram.cpp
	$(CPP) $(FLAGS) -c histogram.cpp

distance.o: distance.cpp
	$(CPP) $(FLAGS) -c distance.cpp

minimization.o: minimization.cpp
	$(CPP) $(FLAGS) -c minimization.cpp

readwritedisplay.o: readwritedisplay.cpp
	$(CPP) $(FLAGS) -c readwritedisplay.cpp

clean:
	del *.o
