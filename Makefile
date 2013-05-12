.PHONY: clean

OBJS    = cells.o list.o table.o pairs.o simulate.o

all: simulate.x

simulate.x: $(OBJS)
	$(CC) $(OBJS) -o $@

clean:
	-rm *.o simulate.x
