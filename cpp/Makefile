# Author: Caroline Moore-Kochlacs, May 2013
# modified from Justin Kinney, Sep 2008

CXX = g++
FLAGS = -Wall -g -std=c++11 -O3
LIBS = 
	
SRCS = simtraces.cc neuron.cc noisenode.cc location.cc direction.cc waveform.cc detector.cc
OBJS = $(SRCS:%.cc=%.o)
DEPSINCLUDE = $(addprefix .deps/,$(SRCS:%.cc=%.d))

default: $(EXECUTABLE)

simtraces: $(OBJS)
	$(CXX) $(FLAGS) $(OBJS) $(LIBS) -o simtraces

%.o: %.cc
	$(CXX) $(FLAGS) -c $< -o $@

clean:
	rm *.o

