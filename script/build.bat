@echo off
REM 0. to the root directory
pushd %~dp0\..

echo Starting build process...

REM Delete build folder
if exist addons\h3-godot\bin rmdir /s /q addons\h3-godot\bin

REM Copy GDExtension configuration file
echo.
echo Copying GDExtension configuration file...
if not exist addons\h3-godot\bin mkdir addons\h3-godot\bin
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
REM Disable docs generation to avoid requiring Doxygen/Graphviz on CI
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=OFF -DENABLE_DOCS=OFF ..
if errorlevel 1 goto error

echo Building h3...
cmake --build . --config Release
if errorlevel 1 goto error

echo Copying h3.dll to addons\h3-godot\bin...
copy /Y bin\Release\h3.dll ..\..\addons\h3-godot\bin\
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

REM remove .exp and .lib file under addons\h3-godot\bin
del addons\h3-godot\bin\*.exp
del addons\h3-godot\bin\*.lib

echo.
echo Build completed successfully!

REM Copy plugin files from export to addons because addons/ is not uploaded
echo.
echo Copying plugin files from export to addons\h3-godot\ ...
if exist export\h3-godot.gd (
	copy /Y export\h3-godot.gd addons\h3-godot\
)
if exist export\plugin.cfg (
	copy /Y export\plugin.cfg addons\h3-godot\
)
if exist export\h3_godot_extension.gdextension (
	copy /Y export\h3_godot_extension.gdextension addons\h3-godot\
)


REM back to the source directory
popd
goto end

:error
echo Build failed!
popd
exit /b 1

:end 