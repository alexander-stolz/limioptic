#!/usr/bin/env python

import urllib as ul
import sys
import os

print "reading file list.."
a = ul.urlopen("http://ams.amstolz.de/update.txt")
b = a.readlines()
a.close()

q = raw_input("press ENTER to start update process. press \"q\" to cancel..  ")

if sys.argv[0].endswith(".py"):
    if (q != "q"):
        for i in range(1, len(b)):
            ul.urlretrieve(b[i].split()[0], b[i].split()[1])
            print "{} updated".format(b[i].split()[1])
        print "\n", b[0]
    try:
        os.remove("liblimioptic-linux.so")
    except:
        pass
else:
    print "precompiled version detected"
    ul.urlretrieve(r"https://github.com/alexander-stolz/limioptic/raw/master/precompiled%20for%20windows/___LIMIOPTIC.exe", "___LIMIOPTIC.exe")
    print "___LIMIOPTIC.exe updated"
    ul.urlretrieve(r"https://github.com/alexander-stolz/limioptic/raw/master/precompiled%20for%20windows/liblimioptic-win.so", "liblimioptic-win.so")
    print "liblimioptic-win updated"

os.System("chmod 777 ___RUN__LIMIOPTIC_2.py")

print "update complete."
raw_input("press enter to exit")
