CPPFLAGS=-std=c++11 -Wall -O2

all: main debug

main: main.cpp TwinGreedyFast.h Graph.h Competitor.h
	g++ main.cpp -o maxcut  ${CPPFLAGS}
debug: main.cpp TwinGreedyFast.h Graph.h Competitor.h
	g++ main.cpp -o maxcut_debug -std=c++11 -Wall -Og -g 