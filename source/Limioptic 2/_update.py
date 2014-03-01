#!/usr/bin/env python

import urllib as ul
import sys

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
	    print "\n",b[0]

	    print "update complete."
	    raw_input("press enter to exit")
else:
	print "precompiled version detected"
	try:
		sys.winver
		print "windows os detected"
		ul.urlretrieve(r"https://github.com/alexander-stolz/limioptic/raw/master/precompiled%20for%20windows/___LIMIOPTIC.exe", "___LIMIOPTIC.exe")
		print "___LIMIOPTIC.exe updated"
		ul.urlretrieve(r"https://github.com/alexander-stolz/limioptic/raw/master/precompiled%20for%20windows/liblimioptic-win.so", "liblimioptic-win.so")
		print "liblimioptic-win updated"
	except:
		print "linux os detected"
		ul.urlretrieve(r"https://github.com/alexander-stolz/limioptic/raw/master/precompiled%20for%20windows/___LIMIOPTIC.exe", "___LIMIOPTIC.exe")
		print "___LIMIOPTIC.exe updated"
		ul.urlretrieve(r"https://github.com/alexander-stolz/limioptic/raw/master/precompiled%20for%20windows/liblimioptic-linux.so", "liblimioptic-linux.so")
		print "liblimioptic-linux updated"