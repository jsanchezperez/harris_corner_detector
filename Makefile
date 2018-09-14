CC=gcc
CXX=g++
CFLAGS=-O3 -Wall -Wextra -fopenmp
LFLAGS=-lpng -ljpeg -ltiff -lm

all: bin obj bin/harris_corner_detector

bin:
	mkdir -p bin

obj:
	mkdir -p obj

obj/iio.o: src/iio.c src/iio.h
	$(CC) -c $< -o $@ -std=c99 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 

obj/gradient.o: src/gradient.cpp src/gradient.h
	$(CXX) -c $< -o $@ -std=c++11 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 

obj/gaussian.o: src/gaussian.cpp src/gaussian.h
	$(CXX) -c $< -o $@ -std=c++11 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 

obj/interpolation.o: src/interpolation.cpp src/interpolation.h
	$(CXX) -c $< -o $@ -std=c++11 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 

obj/zoom.o: src/zoom.cpp src/zoom.h
	$(CXX) -c $< -o $@ -std=c++11 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 

obj/harris.o: src/harris.cpp
	$(CXX) -c $< -o $@ -std=c++11 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 


	
# ------- Main -------
bin/harris_corner_detector: src/main.cpp obj/harris.o obj/iio.o obj/gradient.o obj/gaussian.o obj/interpolation.o obj/zoom.o 
	$(CXX) -std=c++11 -o $@ $^ $(CFLAGS) $(LFLAGS)

clean:
	rm -f obj/*.o bin/harris_corner_detector
	rmdir bin obj
