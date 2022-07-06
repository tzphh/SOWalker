
CC = g++
INCLUDE = -I.
FLAGS = -std=c++11 -lpthread -lortools -fopenmp -Wall -D FASTSKIP -D EXPECT_SCHEDULE 

apps : test/preprocess test/node2vec test/autoregressive test/gen
test/% : test/%.cpp
	@mkdir -p bin/$(@D)
	$(CC) $@.cpp -o bin/$@ $(INCLUDE) $(FLAGS)

clean :
	-rm -rf bin

clear : 
	-rm dataset/*.deg dataset/*.csr dataset/*.beg dataset/*.rat dataset/*.blocks dataset/*.meta

walks:
	-rm dataset/*.walk
