CC = g++
CFLAGS = -Wall -Werror -O3
INCLUDE = -I.


MAIN = main.o

OBJ = CollideStream.o LBField.o

all: $(OBJ) $(MAIN)
	$(CC) -o lbsim -pg $(OBJ) $(MAIN) $(CFLAGS)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $(INCLUDE) -o $*.o $*.cpp

cleanall clean:
	for name in lbsim main.o $(OBJ); do \
	if [ -f $$name ]; then rm $$name; fi; \
	done;
