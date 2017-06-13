cd setup_files
python-2.7.13.msi
setx path "%PATH%;C:\Python27"

REM "setup_files/VCForPython27.exe"

python get-pip.py
python -m pip install numpy-1.13.0+mkl-cp27-cp27m-win32.whl
python -m pip install scipy-0.19.0-cp27-cp27m-win32.whl
python -m pip install PyQt4-4.11.4-cp27-cp27m-win32.whl
python -m pip install matplotlib
python -m pip install pyqtgraph

PAUSE