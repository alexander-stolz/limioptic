Under some Linux systems there is a known bug which causes the vtkRenderWindowInteractor to slow down the program enormously. Plot->PyQtGraph should work until this problem is fixed.

<br>
<b>Installation (Windows):</b>

<b>NEW: You can now just run the precompiled version for Windows. No installation needed!</b><br>
(this might not work on every machine)

This is the web-installer! Full setup: <a href="https://sourceforge.net/projects/limioptic/">https://sourceforge.net/projects/limioptic/</a>

1. Run <b>setup.exe</b>
2. Install the prerequisites in their default directories.
3. <b>Make sure to select the C++ compiler during the MinGW installation!</b>
4. If you have installed LIMIOPTIC in the default directory, you will not be able to update via update.py because of the windows rights management. <b>To fix this, copy the Limioptic 2 folder to another directory</b> (for example to the desktop).
5. Run update.py in the Limioptic 2 directory to update LIMIOPTIC to the newest version.


<br>
<b>Installation (Linux):</b>

<b>NEW: You can now just run the precompiled version for Windows with WINE. No installation needed!</b><br>
(You will have to use the newest (beta-) version of Wine:  <a href="http://www.winehq.org/download/ubuntu">http://www.winehq.org/download/ubuntu</a>)

1. Download the Limioptic 2 folder
2. Limioptic imports the following libraries. Make sure you have these:
  - python-qt4
  - (python-vtk)
  - <a href="http://www.pyqtgraph.org/">PyQtGraph</a> (depends on numpy and scipy)
  - scipy
  - pylab



<br>
<b>Licence</b>

The program Limioptic 2 maintained by Alexander Stolz is freely available and distributable. However, if you use it for some work whose results are made public, then you have to reference it properly.
