#!/usr/bin/env python

import urllib.request as ul
import sys
import os

print("reading file list..")
a = ul.urlopen("https://raw.githubusercontent.com/alexander-stolz/limioptic/master/LIMIOPTIC/update.txt").read()
b = "".join(map(chr, a)).split("\n")
print(b)

q = input("press ENTER to start update process. press \"q\" to cancel..  ")

if sys.argv[0].endswith(".py"):
    if (q != "q"):
        for i in range(1, len(b)):
            try:
                ul.urlretrieve(b[i].split()[0], b[i].split()[1])
                print("{} updated".format(b[i].split()[1]))
            except:
                try:
                    print(b[i].split()[1], "not found")
                except:
                    pass
        print("\n", b[0])
    try:
        os.remove("liblimioptic-linux.so")
    except:
        pass

    try:
        os.remove("___RUN__LIMIOPTIC_2.py")
        os.remove("___RUN__LIMIOPTIC_2.bat")
    except:
        pass


else:
    print("precompiled version detected")
    ul.urlretrieve(r"https://github.com/alexander-stolz/limioptic/raw/master/precompiled%20for%20windows/___LIMIOPTIC.exe", "___LIMIOPTIC.exe")
    print("___LIMIOPTIC.exe updated")
    ul.urlretrieve(r"https://github.com/alexander-stolz/limioptic/raw/master/precompiled%20for%20windows/liblimioptic-win.so", "liblimioptic-win.so")
    print("liblimioptic-win updated")
    ul.urlretrieve(r"https://github.com/alexander-stolz/limioptic/raw/master/precompiled%20for%20windows/_Calculator.exe", "_Calculator.exe")
    print("_Calculator.exe updated")

print("update complete.")
input("press enter to exit")
