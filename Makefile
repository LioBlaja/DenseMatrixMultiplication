CC = gcc
OBJ = func

all:
	# $(CC) -fopenmp -O3 -DDEBUG -Wall -O3 main.c ./permutedMultiplications/serialPMA.c ./permutedMultiplications/parallelPMA.c -o $(OBJ)
	$(CC) -fopenmp -O3 -Wall  main.c ./permutedMultiplications/serialPMA.c ./permutedMultiplications/parallelPMA.c  -o $(OBJ)
exec:
	# ./$(OBJ) > file.txt
	./$(OBJ)
clean:
	rm -rf $(OBJ)
