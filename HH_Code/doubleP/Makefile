CC=g++
SRCS=src/*.cpp
FLAGS=--std=c++11 -O2
INCLUDES=-I./includes

ifndef DT
	DT=0.01
	PREREQ=$(SRCS) includes/*
else
	PREREQ=FORCE
endif

all: runner

runner: $(PREREQ)
	$(CC) $(FLAGS) $(INCLUDES) $(SRCS) -DDT=$(DT) -o runner

.PHONY: run plot clean FORCE

run: runner
	mkdir -p results
	./runner -case 21 -nBurst 200 -pExcN 1.0 -vInh 70 > ./results/dt$(subst 0.,,$(DT)).out
	rm ./results/*.txt

plot:
	$(MAKE) -C ../ plot

clean:
	rm runner

FORCE:
