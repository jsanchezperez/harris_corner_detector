CC=gcc
CXX=g++
CFLAGS=-Wall -Wextra -O3 -Werror -fopenmp
LFLAGS=-lpng -ljpeg -ltiff -lm

all: bin obj bin/harris_corner_detector

bin:
	mkdir -p bin

obj:
	mkdir -p obj

obj/iio.o: src/iio.c src/iio.h
	$(CC) -c $< -o $@ -std=c99 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 

obj/harris.o: src/harris.cpp
	$(CXX) -c $< -o $@ -std=c++11 $(CFLAGS) -Wno-unused -pedantic -DNDEBUG -D_GNU_SOURCE 


	
# ------- Main -------
bin/harris_corner_detector: src/main.cpp obj/harris.o obj/iio.o
	$(CXX) -std=c++11 -o $@ $^ $(CFLAGS) $(LFLAGS)

	
clean:
	rm -f obj/*.o bin/harris_corner_detector
	rmdir bin obj
