@echo off
REM 0. to the root directory
pushd %~dp0\..

echo Starting cleanup process...

REM Clean godot-cpp
echo.
echo Cleaning godot-cpp...
scons -C godot-cpp --clean
if errorlevel 1 goto error

REM Clean h3 library
echo.
echo Cleaning h3 library...
if exist h3\build (
    echo Removing h3/build directory...
    rmdir /s /q h3\build
)

REM Clean GDExtension
echo.
echo Cleaning GDExtension...
if exist build (
    echo Removing build directory...
    rmdir /s /q build
)

REM Clean obj files in src directory
echo.
echo Cleaning obj files in src directory...
for /r src %%F in (*.obj) do (
    echo Removing %%F
    del /f /q "%%F"
)

REM Clean gen directory
echo.
echo Cleaning gen directory...
if exist src\gen (
    echo Removing src/gen directory...
    rmdir /s /q src\gen
)

REM Clean SCons database file
if exist .sconsign.dblite (
    echo Removing SCons database file...
    del /f /q .sconsign.dblite
)

echo.
echo Cleanup completed successfully!

REM back to the source directory
popd
goto end

:error
echo Cleanup failed!
popd
exit /b 1

:end 