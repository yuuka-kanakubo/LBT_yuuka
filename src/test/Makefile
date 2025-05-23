# Makefile for LBT energy-momentum conservation test

CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall
INCLUDES = -I./include
SRCS = test_LBTEMconservation.cpp ../LBTcl.cpp ../LBTConfig.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = test_LBT_energyconservation

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

