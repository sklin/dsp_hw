.PHONY: all clean

CFLAGS+=
CXXFLAGS+=-std=c++11 -g
LDFLAGS+=-lm     # link to math library

TARGET=test train

all: $(TARGET)
# type make/make all to compile test_hmm

test: test_hmm.c hmm.h
	g++ test_hmm.c -o test

train: train_hmm.cpp hmm.h
	g++ train_hmm.cpp -o train -std=c++11


clean:
	$(RM) $(TARGET)   # type make clean to remove the compiled file
