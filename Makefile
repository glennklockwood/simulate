.PHONY: clean

CC      = c99
CFLAGS  = -g -O0

OBJS    = cells.o list.o table.o pairs.o simulate.o move.o

all: simulate.x

simulate.x: $(OBJS)
	$(CC) $(OBJS) -o $@ -lm

clean:
	-rm *.o simulate.x
