.PHONY: all

all: base optimizer1 optimizer2 optimizer3 optimizer4

base: base.cpp
	mpic++ -std=c++11 ./base.cpp -o base

optimizer1: optimizer1.cpp
	mpic++ -std=c++11 ./optimizer1.cpp -o optimizer1

optimizer2: optimizer2.cpp
	mpic++ -std=c++11 ./optimizer2.cpp -o optimizer2

optimizer3: optimizer3.cpp
	mpic++ -std=c++11 ./optimizer3.cpp -o optimizer3

optimizer4: optimizer4.cpp
	mpic++ -std=c++11 ./optimizer4.cpp -o optimizer4

clean:
	rm base optimizer1 optimizer2 optimizer3 optimizer4