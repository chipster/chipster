@echo off
setlocal

rem Copyright (c) 1999, 2010 Tanuki Software, Ltd.
rem http://www.tanukisoftware.com
rem All rights reserved.
rem
rem This software is the proprietary information of Tanuki Software.
rem You shall use it only in accordance with the terms of the
rem license agreement you entered into with Tanuki Software.
rem http://wrapper.tanukisoftware.com/doc/english/licenseOverview.html
rem
rem Java Service Wrapper command based script.

rem -----------------------------------------------------------------------------
rem These settings can be modified to fit the needs of your application
rem Optimized for use with version 3.5.7 of the Wrapper.

rem The base name for the Wrapper binary.
set _WRAPPER_BASE=wrapper

rem The name and location of the Wrapper configuration file.
rem  (Do not remove quotes.)
set _WRAPPER_CONF="../conf/wrapper.conf"

rem _FIXED_COMMAND tells the script to use a hard coded action rather than
rem  expecting the first parameter of the command line to be the command.
rem  By default the command will will be expected to be the first parameter.
rem set _FIXED_COMMAND=console

rem _PASS_THROUGH tells the script to pass all arguments through to the JVM
rem  as is.  If _FIXED_COMMAND is specified then all arguments will be passed.
rem  If not set then all arguments starting with the second will be passed.
set _PASS_THROUGH=true

rem Do not modify anything beyond this point
rem -----------------------------------------------------------------------------

if "%OS%"=="Windows_NT" goto nt
echo This script only works with NT-based versions of Windows.
goto :eof

:nt
rem Find the application home.
rem %~dp0 is location of current script under NT
set _REALPATH=%~dp0

rem
rem Decide on the specific Wrapper binary to use (See delta-pack)
rem
if "%PROCESSOR_ARCHITECTURE%"=="AMD64" goto amd64
if "%PROCESSOR_ARCHITECTURE%"=="IA64" goto ia64
set _WRAPPER_L_EXE=%_REALPATH%%_WRAPPER_BASE%-windows-x86-32.exe
goto search
:amd64
set _WRAPPER_L_EXE=%_REALPATH%%_WRAPPER_BASE%-windows-x86-64.exe
goto search
:ia64
set _WRAPPER_L_EXE=%_REALPATH%%_WRAPPER_BASE%-windows-ia-64.exe
goto search
:search
set _WRAPPER_EXE=%_WRAPPER_L_EXE%
if exist "%_WRAPPER_EXE%" goto validate
set _WRAPPER_EXE=%_REALPATH%%_WRAPPER_BASE%.exe
if exist "%_WRAPPER_EXE%" goto validate
echo Unable to locate a Wrapper executable using any of the following names:
echo %_WRAPPER_L_EXE%
echo %_WRAPPER_EXE%
pause
goto :eof

:validate
if not [%_FIXED_COMMAND%]==[] goto defaultaction
rem There should be a command on the command line.  Look for it.
set _COMMAND=
for /F %%v in ('echo %1^|findstr "^console$ ^start$ ^pause$ ^resume$ ^stop$ ^restart$ ^install$ ^remove$"') do call :exec set _COMMAND=%%v

if [%_COMMAND%]==[] (
    set _COMMAND=%1
    goto showusage
) else (
    rem Got a command
    shift
)
goto havecommand
    
:defaultaction	
rem Specified a default action.
set _COMMAND=%_FIXED_COMMAND%
goto havecommand

:havecommand
if [%_PASS_THROUGH%]==[] goto callcommand
rem Collect an parameters
:parameters
set _PARAMETERS=%_PARAMETERS% %1
shift
if not [%1]==[] goto parameters

:callcommand
rem
rem Run the application.
rem At runtime, the current directory will be that of wrapper.exe
rem
call :%_COMMAND%
if errorlevel 1 goto handleerror
goto :eof

:handleerror
if [%_MATCHED]==[] goto showusage
pause
goto :eof

:showusage
rem A command was not specified, or it was now known.
if not [%_COMMAND%]==[] (
    echo Unknown command: %_COMMAND%
    echo.
)
if [%_PASS_THROUGH%]==[] (
    echo Usage: %0 [ console : start : pause : resume : stop : restart : install : update : remove ]
) else (
    echo Usage: %0 [ console {JavaAppArgs} : start : pause : resume : stop : restart : install {JavaAppArgs} : update {JavaAppArgs} : remove ]
)
pause
goto :eof

:console
set _MATCHED=true
if [%_PASS_THROUGH%]==[] (
    "%_WRAPPER_EXE%" -c %_WRAPPER_CONF%
) else (
    "%_WRAPPER_EXE%" -c %_WRAPPER_CONF% -- %_PARAMETERS%
)
goto :eof

:start
set _MATCHED=true
"%_WRAPPER_EXE%" -t %_WRAPPER_CONF%
goto :eof

:pause
set _MATCHED=true
"%_WRAPPER_EXE%" -a %_WRAPPER_CONF%
goto :eof

:resume
set _MATCHED=true
"%_WRAPPER_EXE%" -e %_WRAPPER_CONF%
goto :eof

:stop
set _MATCHED=true
"%_WRAPPER_EXE%" -p %_WRAPPER_CONF%
goto :eof

:install
set _MATCHED=true
if [%_PASS_THROUGH%]==[] (
    "%_WRAPPER_EXE%" -i %_WRAPPER_CONF%
) else (
    "%_WRAPPER_EXE%" -i %_WRAPPER_CONF% -- %_PARAMETERS%
)
goto :eof

:remove
set _MATCHED=true
"%_WRAPPER_EXE%" -r %_WRAPPER_CONF%
goto :eof

:restart
set _MATCHED=true
call :stop
call :start
goto :eof

:exec
%*
goto :eof
