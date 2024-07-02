CC = g++
CFLAGS = -O3 -std=c++11

btl_b: main.o hung.o
	$(CC) $(CFLAGS) -o btl main.o hung.o

hung.o: Hungarian.cpp Hungarian.h
	$(CC) -c Hungarian.cpp -o hung.o
	
main.o: btl.cpp Hungarian.h
	$(CC) $(CFLAGS) -c btl.cpp -o main.o

clean:
	-rm main.o hung.o

