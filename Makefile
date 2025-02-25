CC = /opt/homebrew/bin/gcc-14  # Update with the correct path
CFLAGS = -O2 -fopenmp
LDFLAGS = -lm -fopenmp
SRC = solver.c Lab3IO.c
OBJ = $(SRC:.c=.o)
EXEC = main

DATA = datagen.c
DATA_EXEC = datagen

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJ) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXEC)

data: 
	$(CC) -o $(DATA_EXEC) $(DATA)
	./datagen

run: $(EXEC)
	./$(EXEC)