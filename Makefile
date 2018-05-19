CC=gcc
CFLAGS=
OBJECTS=molecular_simulation.c particles.o files.o

nn_cc: $(OBJECTS)
	$(CC) -o ./body3 $(OBJECTS) $(CFLAGS)

%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

clean:
	rm -f *.o body3
