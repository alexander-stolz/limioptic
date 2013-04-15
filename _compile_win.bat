@echo This will only work if MinGW is installed!
@PAUSE
@echo.
@echo compiling climioptic.cpp
@echo ------------------------
@g++ -Wall -fPIC -c climioptic.cpp
@echo.
@echo.
@echo compiling limioptic.cpp
@echo ------------------------
@g++ -Wall -fPIC -c limioptic.cpp
@echo.
@echo.
@echo building library
@echo ------------------------
@g++ -fPIC -shared -Wl,-soname,liblimioptic.so -o liblimioptic.so limioptic.o climioptic.o
@echo completed.