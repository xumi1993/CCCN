OBJ = sacio.o stack.o
stack : $(OBJ)
	cc -o ../bin/stack $(OBJ)
$(OBJ) : sacio.h
sacio.o stack.o : sacio.h

clean:
	rm -f $(OBJ)
