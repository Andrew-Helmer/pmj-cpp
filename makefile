CXXFLAGS = -Wall -std=c++17

all: generate_samples

release: CXXFLAGS += -g3
release: generate_samples

DEBUG: debug

debug: CXXFLAGS += -DDEBUG -g
debug: generate_samples

generate_samples: generate_samples.cc pj.cc pmj.cc pmj02.cc util.cc
	g++ $(CXXFLAGS) -o generate_samples generate_samples.cc pj.cc pmj.cc pmj02.cc util.cc -I

clean:
	rm -Rf generate_samples generate_samples.dSYM