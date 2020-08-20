powershell -Command "Invoke-WebRequest https://www.python.org/ftp/python/3.8.5/python-3.8.5.exe -OutFile python-3.8.5.exe"
python-3.8.5.exe PrependPath=1
del python-3.8.5.exe