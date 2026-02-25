CXX ?= g++-11
CC  ?= gcc-11

CFLAGS   := -O3 -ggdb -g3
CXXFLAGS := -O3 -ggdb -g3

ifneq ($(DBG),)
	CFLAGS   += -DGFA_ED_DBG=4
	CXXFLAGS += -DGFA_ED_DBG=4
endif
ifneq ($(GDB),)
	CFLAGS   += -O0 -ggdb
	CXXFLAGS += -O0 -ggdb
endif
ifneq ($(ASAN),)
	CFLAGS   += -fsanitize=address
	CXXFLAGS += -fsanitize=address
	LDFLAGS  += -fsanitize=address
endif

all: gwfa

gwfa.o: gwfa.c gwfa.h ksort.h kvec.h
	$(CC) $(CFLAGS) -c -o $@ $<

main.o: main.cpp gwfa.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

gwfa: main.o gwfa.o
	$(CXX) $(LDFLAGS) -o $@ $^

clean:
	rm -f *.o gwfa
