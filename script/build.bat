@echo off
REM 0. to the root directory
pushd %~dp0\..

echo Starting build process...

REM Delete build folder
if exist build rmdir /s /q build

REM Copy GDExtension configuration file
echo.
echo Copying GDExtension configuration file...
if not exist build\bin mkdir build\bin
copy /Y h3_godot_extension.gdextension build\bin\
copy /Y h3_godot_extension.gdextension.uid build\bin\
if errorlevel 1 goto error

REM Update and build godot-cpp
echo.
echo Building godot-cpp...
git submodule update --init --recursive

echo Building godot-cpp release version...
scons -C godot-cpp target=template_release use_static_cpp=yes
if errorlevel 1 goto error

REM Build h3 library
echo.
echo Building h3 library...

pushd h3
if not exist build mkdir build
cd build

echo Configuring CMake...
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=OFF ..
if errorlevel 1 goto error

echo Building h3...
cmake --build . --config Release
if errorlevel 1 goto error

echo Copying h3.dll to build/bin...
copy /Y bin\Release\h3.dll ..\..\build\bin\
if errorlevel 1 goto error

popd

REM Build GDExtension
echo.
echo Building GDExtension...

echo Building GDExtension debug version...
scons target=template_debug platform=windows
@REM scons target=template_debug platform=linux
@REM scons target=template_debug platform=macos
@REM scons target=template_debug platform=ios
@REM scons target=template_debug platform=android
@REM scons target=template_debug platform=web
if errorlevel 1 goto error

echo Building GDExtension release version...
scons target=template_release platform=windows
@REM scons target=template_release platform=linux
@REM scons target=template_release platform=macos
@REM scons target=template_release platform=ios
@REM scons target=template_release platform=android
@REM scons target=template_release platform=web
if errorlevel 1 goto error

echo.
echo Build completed successfully!
echo Debug version: build/bin/libh3_extension.windows.template_debug.x86_64.dll
echo Release version: build/bin/libh3_extension.windows.template_release.x86_64.dll

REM back to the source directory
popd
goto end

:error
echo Build failed!
popd
exit /b 1

:end 