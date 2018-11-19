COMP = g++ -std=c++11 -c -g
LINK = g++ -std=c++11 -lpthread

all: main
clean:
	rm -f *.o rm main
main: main.o
	$(LINK) main.o -o main
main.o: main.cpp first_lower_bound.h mathFunctions.h random_points.h reading.h
	$(COMP) main.cpp -o main.o;
