#!/bin/make

TARGET = ./main$(EXEEXT)
#SRCS = main.cpp vector3.cpp screen.cpp color.cpp light.cpp photon.cpp sphere.cpp polygon3.cpp scene.cpp aabb3.cpp aabb3n.cpp
#SRCS = ulist.c hamming.c main.c
SRCS = unit.c hamming.c main.c
CXX = gcc
CXXFLAGS = -g -Wall # -g

#OPENCVINC = `pkg-config --cflags opencv`
#OPENCVLIB = `pkg-config --libs opencv`

#LDFLAGS  = -fopenmp
LDFLAGS  = 

OBJS = $(SRCS:.c=.o)

.PHONY: all clean
.SUFFIXES: .c .o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS) $(OPENCVLIB)

.c.o:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OPENCVINC) -c $<

depend:
	$(CXX) -MM $(INCLUDE) $(CXXFLAGS) $(SRCS) > dependencies

clean:
	rm -rf $(OBJS) $(TARGET)


