#!/bin/sh

echo "g++ -Wall -fPIC -O -c climioptic.cpp"
g++ -std=c++0x -Wall -fPIC -O -c climioptic.cpp

echo "g++ -Wall -fPIC -O -c limioptic.cpp"
g++ -Wall -fPIC -O -c limioptic.cpp

echo "g++ -fPIC -shared -Wl,-soname,liblimioptic.so -O -o liblimioptic.so limioptic.o climioptic.o"
g++ -fPIC -shared -Wl,-soname,liblimioptic-linux.so -O -o liblimioptic-linux.so limioptic.o climioptic.o

