CC = g++
CFLAGS = -g -std=c++11
LFLAGS = 
TARGET = sph

all: $(TARGET)

$(TARGET): sph.o main.o
	g++ $(LFLAGS) sph.o main.o -o sph

sph.o: sph.h sph.C
	g++ $(CFLAGS) -c sph.C

main.o: sph.h main.C
	g++ $(CFLAGS) -c main.C

clean:
	$(RM) -f $(TARGET) *.o *.out
