#!/bin/sh

echo "g++ -Wall -fPIC -O -c climioptic.cpp"
g++ -Wall -fPIC -O -c climioptic.cpp

echo "g++ -Wall -fPIC -O -c limioptic.cpp"
g++ -Wall -fPIC -O -c limioptic.cpp

echo "g++ -fPIC -shared -Wl,-soname,liblimioptic.so -O -o liblimioptic.so limioptic.o climioptic.o"
g++ -fPIC -shared -Wl,-soname,liblimioptic.so -O -o liblimioptic.so limioptic.o climioptic.o

