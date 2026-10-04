@echo off
rem Windows shim for the fake genome-cli: test/bin/genome is a sh script.
rem The tests set GENETICS_FAKE_GENOME_SH to sh.exe (Git for Windows, MSYS2).
"%GENETICS_FAKE_GENOME_SH%" "%~dp0genome" %*
exit /b %ERRORLEVEL%
