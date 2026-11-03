MODULE_DIRECTORY ?= ../cpp_modules

OBJS := build/src/test.o 
DEPS := $(OBJS:.o=.d)

all: LDFLAGS= -fuse-ld=gold -L../ -lz -ldl -pthread
all: build/src/test.o
	@mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $(wildcard build/src/*.o) -o test 

build/src/test.o: CXXFLAGS= -O3 -pthread -std=c++23 -I../include -I$(MODULE_DIRECTORY)/expected/include -I$(MODULE_DIRECTORY)/graph-v3/include -I$(MODULE_DIRECTORY)/doctest  -MMD
build/src/test.o: src/test.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c -o build/src/test.o src/test.cpp

clean:
	cd test/build && $(MAKE) clean

-include $(DEPS)
