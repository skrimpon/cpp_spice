#
# Skrimponis Panagiotis

CC = g++ 
LD_FLAGS =
CC_FLAGS = -w -O3 -std=c++11 -Wall -lm

EXEC = cppspice
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(patsubst src/%,build/%,$(SOURCES:.cpp=.o))
INCLUDE = -I include

$(EXEC): $(OBJECTS)
	$(CC) $(LD_FLAGS) $(OBJECTS) -o $(EXEC) /usr/local/lib/libcxsparse.a

build/%.o: src/%.cpp
	$(CC) -c $(INCLUDE) $(CC_FLAGS) $< -o $@

clean:
	rm -rf $(EXEC) $(OBJECTS) *.txt

