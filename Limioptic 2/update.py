#!/usr/bin/env python

import urllib as ul

print "reading file list.."
a = ul.urlopen("http://ams.amstolz.de/update.txt")
b = a.readlines()
a.close()

q = raw_input("press ENTER to start update process. press \"q\" to cancel..  ")

if (q != "q"):
    for i in range(1,len(b)):
        ul.urlretrieve(b[i].split()[0],b[i].split()[1])
        print "{} updated".format(b[i].split()[1])
    print "\n",b[0]

    print "update complete."
    raw_input("press enter to exit")
