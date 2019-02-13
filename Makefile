CFLAGS = --std=gnu++11 -Wall -Wextra -Ofast
CC = g++ $(CFLAGS)

unaryAutomata: unaryAutomata.cpp unaryAutomata.hpp trees.hpp constants.hpp
	$(CC) -o unaryAutomata unaryAutomata.cpp

clean:
	rm -f *.o

distclean: clean
	rm -f unaryAutomata
