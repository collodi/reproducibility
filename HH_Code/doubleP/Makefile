CC=g++
SRCS=src/*.cpp
FLAGS=--std=c++11
INCLUDES=-I./includes
DT=0.01

all: runner

runner: $(SRCS)
	$(CC) $(FLAGS) $(INCLUDES) $(SRCS) -Ddt=$(DT) -o runner

.PHONY: run clean

run: runner
	mkdir -p results
	./runner -case 21 -nBurst 200 -pExcN 1.0 -vInh 70

clean:
	rm runner
