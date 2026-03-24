#!/usr/bin/env bash
set -euo pipefail

# Switch to repo root (script is in script/)
cd "$(dirname "$0")/.."

echo "Starting build process..."

# Remove old build output
echo "Removing old addons/h3-godot/bin..."
rm -rf addons/h3-godot/bin

echo "Creating addons/h3-godot/bin..."
mkdir -p addons/h3-godot/bin

echo
echo "Updating git submodules..."
git submodule update --init --recursive

echo
echo "Building godot-cpp (release)..."
scons -C godot-cpp target=template_release use_static_cpp=yes

echo
echo "Building h3 library..."
pushd h3 > /dev/null
mkdir -p build
cd build

echo "Configuring CMake..."
# Disable building docs in CI/automation to avoid requiring Doxygen/Graphviz on runners
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=OFF -DENABLE_DOCS=OFF ..

echo "Building h3..."
cmake --build . --config Release

# Try to locate the built library (support common linux names and fallbacks)
LIB=$(find . -type f \( -name "libh3*.so*" -o -name "h3*.so*" -o -name "libh3*.a" -o -name "h3*.a" \) | head -n1 || true)
if [ -z "$LIB" ]; then
  if [ -f bin/Release/h3.dll ]; then
    LIB=bin/Release/h3.dll
  fi
fi
if [ -z "$LIB" ]; then
  echo "Built library not found in h3/build. Listing contents:" >&2
  ls -la || true
  popd > /dev/null
  exit 1
fi

echo "Copying $LIB to addons/h3-godot/bin..."
cp -v "$LIB" ../../addons/h3-godot/bin/

# create unversioned soname symlink in addons bin so .gdextension dependency matches
BASENAME=$(basename "$LIB")
if [[ "$BASENAME" == libh3.so.* || "$BASENAME" == libh3.so ]]; then
  (cd ../../addons/h3-godot/bin && ln -sf "$BASENAME" libh3.so) || true
fi

popd > /dev/null

echo
# Determine architecture for naming (default to x86_64)
ARCH=$(uname -m)
case "$ARCH" in
  x86_64) ARCH_NAME="x86_64" ;;
  i386|i686) ARCH_NAME="x86_32" ;;
  aarch64) ARCH_NAME="arm64" ;;
  armv7l|armv6l) ARCH_NAME="arm32" ;;
  *) ARCH_NAME="$ARCH" ;;
esac

echo "Building GDExtension (debug)..."
scons target=template_debug platform=linux

# After debug build, copy/rename extension to match .gdextension naming
DEBUG_OUT=addons/h3-godot/bin/libh3_extension.linux.template_debug.${ARCH_NAME}.so
if [ -f addons/h3-godot/bin/libh3_extension.so ]; then
  cp -v addons/h3-godot/bin/libh3_extension.so "$DEBUG_OUT"
fi

echo "Building GDExtension (release)..."
scons target=template_release platform=linux

# After release build, copy/rename extension to match .gdextension naming
RELEASE_OUT=addons/h3-godot/bin/libh3_extension.linux.template_release.${ARCH_NAME}.so
if [ -f addons/h3-godot/bin/libh3_extension.so ]; then
  cp -v addons/h3-godot/bin/libh3_extension.so "$RELEASE_OUT"
fi

echo
echo "Cleaning Windows artifacts if present..."
rm -f addons/h3-godot/bin/*.exp addons/h3-godot/bin/*.lib || true

echo
echo "Build completed successfully!"

echo "Copying plugin files to addons/h3-godot/..."
cp -v export/h3-godot.gd addons/h3-godot/ || true
cp -v export/plugin.cfg addons/h3-godot/ || true
cp -v export/h3_godot_extension.gdextension addons/h3-godot/ || true

echo
echo "Note: make executable with: chmod +x script/build.sh"
