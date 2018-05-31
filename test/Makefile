CC=g++

all: global stream nestfunc

global: global.cpp
	$(CC) global.cpp -o global

stream: stream.cpp
	$(CC) -std=c++11 stream.cpp -o stream

nestfunc: nestfunc.cpp
	$(CC) nestfunc.cpp -o nestfunc

.PHONY: clean
clean:
	rm global nestfunc stream
