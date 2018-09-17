CC=gcc
CXX=g++
CFLAGS=-O3 -Wall -Wextra -fopenmp 
LFLAGS=-lpng -ljpeg -ltiff -lm 
CVFLAGS=$(shell pkg-config --cflags --libs opencv)


all: bin bin/test

bin:
	mkdir -p bin
	mkdir -p tmp
	mkdir -p obj

obj/iio.o: ../src/iio.c ../src/iio.h
	$(CC) -c $< -o $@ -std=c99 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 

obj/gradient.o: ../src/gradient.cpp ../src/gradient.h
	$(CXX) -c $< -o $@ -std=c++11 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 

obj/gaussian.o: ../src/gaussian.cpp ../src/gaussian.h
	$(CXX) -c $< -o $@ -std=c++11 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 

obj/interpolation.o: ../src/interpolation.cpp ../src/interpolation.h
	$(CXX) -c $< -o $@ -std=c++11 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 

obj/harris.o: ../src/harris.cpp
	$(CXX) -c $< -o $@ -std=c++11 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 

obj/bicubic_interpolation.o: bicubic_interpolation.cpp
	$(CXX) -c $< -o $@ -std=c++11 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 

obj/zoom.o: ../src/zoom.cpp ../src/zoom.h
	$(CXX) -c $< -o $@ -std=c++11 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 

obj/harris_opencv.o: harris_opencv.cpp
	$(CXX) -c $< -o $@ -std=c++11 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 


	
# ------- test -------
bin/test: test.cpp obj/bicubic_interpolation.o obj/harris.o obj/iio.o obj/gradient.o obj/gaussian.o obj/interpolation.o  obj/zoom.o  obj/harris_opencv.o
	$(CXX) -std=c++11 -o $@ $^ $(CFLAGS) $(LFLAGS) $(CVFLAGS)

	
clean:
	rm -f obj/*.o bin/test tmp/*
	rm -R bin obj tmp
