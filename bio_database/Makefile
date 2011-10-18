CC = g++

ifeq ($(shell uname -s), Darwin)
	CFLAGS = -I../berkeley_db/build_osx
	LNFLAGS = -L../berkeley_db/build_osx -lpthread -ldb
else
	CFLAGS = -I../berkeley_db/build_unix
	LNFLAGS = -L../berkeley_db/build_unix -lpthread -ldb
endif


# find all the source files 
SRC = biodb.c
#$(shell ls *.c)
#SRC += $(shell ls */*.cpp)

OBJ=$(SRC:.c=.o)

%.o: %.cpp
	$(CC) -c $< $(CFLAGS)

all: $(OBJ)
	$(CC) -o test $(OBJ) $(LNFLAGS) 

clean:
	rm *.o test