cd setup_files
@echo.
@echo Please install python to the default directory C:\Python27. Otherwise you have to change the windows path environment variable accordingly.
@echo.
@PAUSE
python-2.7.13.msi
setx path "%PATH%;C:\Python27"

REM "setup_files/VCForPython27.exe"

