    rcpgenerator.exe: rcpgenerator.o distance.o minimization.o histogram.o readwritedisplay.o twister.o
	g++ rcpgenerator.o distance.o minimization.o histogram.o readwritedisplay.o twister.o -o rcpgenerator.exe

    rcpgenerator.o: rcpgenerator.cpp
	g++ -c rcpgenerator.cpp

    twister.o: twister.cpp
	g++ -c twister.cpp

    histogram.o: histogram.cpp
	g++ -c histogram.cpp

    distance.o: distance.cpp
	g++ -c distance.cpp

    minimization.o: minimization.cpp
	g++ -c minimization.cpp

    readwritedisplay.o: readwritedisplay.cpp
	g++ -c readwritedisplay.cpp

