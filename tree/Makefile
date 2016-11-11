include Makefile.in

.PHONY: all exe clean

CC=gcc
TARGET=sph

exe: $(TARGET)
all: exe

# =======
$(TARGET): sph.o params.o state.o interact.o leapfrog.o io_bin.o timing.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

sph.o: sph.c params.h state.h interact.h leapfrog.h io.h timing.h

params.o: params.c params.h
state.o: state.c state.h
interact.o: interact.c interact.h state.h params.h timing.h
leapfrog.o: leapfrog.c leapfrog.h state.h params.h
io_txt.o: io_txt.c io.h
io_bin.o: io_bin.c io.h

%.o: %.c
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $<

io_bin.o: io_bin.c
	$(CC) -c $(CFLAGS) $<

# =======
view: 
	java -jar ../jbouncy/Bouncy.jar run.out

# =======
clean:
	rm -f $(TARGET) *~ *.o run.out
