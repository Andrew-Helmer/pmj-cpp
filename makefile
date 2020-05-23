CXXFLAGS = -Wall -std=c++17

all: pmj

release: CXXFLAGS += -g3
release: pmj

DEBUG: debug

debug: CXXFLAGS += -DDEBUG -g
debug: pmj

pmj: main.cc pj.cc pmj.cc util.cc
	g++ $(CXXFLAGS) -o pmj main.cc pj.cc pmj.cc util.cc -I

clean:
	rm -Rf pmj pmj.dSYM