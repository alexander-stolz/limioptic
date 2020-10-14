@echo This will only work if MinGW is installed!
@PAUSE
@echo.
@echo compiling climioptic.cpp
@echo ------------------------
@g++ -std=c++0x -Wall -fPIC -O -c climioptic.cpp 
@echo.
@echo.
@echo compiling limioptic.cpp
@echo ------------------------
@g++ -Wall -fPIC -O -c limioptic.cpp
@echo.
@echo.
@echo building library
@echo ------------------------
@g++ -fPIC -shared -Wl,-soname,liblimioptic-win.so -O -o liblimioptic-win.so limioptic.o climioptic.o -static-libstdc++
@echo completed.
@echo.
@PAUSE