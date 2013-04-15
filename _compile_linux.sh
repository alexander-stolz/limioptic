#!/bin/sh

echo "g++ -Wall -fPIC -c climioptic.cpp"
g++ -Wall -fPIC -c climioptic.cpp

echo "g++ -Wall -fPIC -c limioptic.cpp"
g++ -Wall -fPIC -c limioptic.cpp

echo "g++ -fPIC -shared -Wl,-soname,liblimioptic.so -o liblimioptic.so limioptic.o climioptic.o"
g++ -fPIC -shared -Wl,-soname,liblimioptic.so -o liblimioptic.so limioptic.o climioptic.o

