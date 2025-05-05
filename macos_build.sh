#!/bin/bash
# Comprehensive macOS build and setup script for iterateKT
# This script handles all macOS-specific fixes in one place

# Exit on error
set -e

# Make script directory the working directory
cd "$(dirname "$0")"

# Set project root
export ITERATEKT="$(pwd)"
echo "Setting ITERATEKT=$ITERATEKT"

# Check for required tools
if ! command -v cmake &>/dev/null; then
    echo "ERROR: cmake not found. Please install it with: brew install cmake"
    exit 1
fi

if ! command -v root-config &>/dev/null; then
    echo "ERROR: ROOT not found. Please install it with: brew install root"
    exit 1
fi

if ! command -v brew &>/dev/null; then
    echo "ERROR: Homebrew not found. Please install it first."
    exit 1
fi

# Check for ninja and install if needed
if ! command -v ninja &>/dev/null; then
    echo "Installing ninja build system..."
    brew install ninja
fi

# Make sure we have Boost
if ! brew list --formula | grep -q boost; then
    echo "Installing Boost..."
    brew install boost
fi

# Create build directory
echo "Creating clean build directory..."
rm -rf build
mkdir -p build
cd build

# Configure with ninja
echo "Configuring with CMake and Ninja..."
cmake -G Ninja \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_FLAGS="-fPIC -O2" \
    -DCMAKE_INSTALL_BINDIR=${ITERATEKT}/bin \
    -DCMAKE_INSTALL_LIBDIR=${ITERATEKT}/lib \
    -DCMAKE_MACOSX_RPATH=ON \
    ..

# Build with ninja
echo "Building with ninja..."
ninja

# Create directories if they don't exist
mkdir -p ${ITERATEKT}/lib
mkdir -p ${ITERATEKT}/bin

# Install to bin and lib directories
echo "Installing..."
ninja install

# Create symbolic link for .so to .dylib (for ROOT compatibility)
cd ${ITERATEKT}/lib
if [ -f "libITERATEKT.dylib" ] && [ ! -f "libITERATEKT.so" ]; then
    echo "Creating symbolic link from .dylib to .so..."
    ln -sf libITERATEKT.dylib libITERATEKT.so
fi

# Return to project root
cd ${ITERATEKT}

# Create run script with proper environment variables
cat > ${ITERATEKT}/run_iterateKT.sh << 'EOF'
#!/bin/bash
# Script to run iterateKT with proper environment variables

# Get directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

# Set environment variables
export ITERATEKT="${SCRIPT_DIR}"
export DYLD_LIBRARY_PATH="${SCRIPT_DIR}/lib:${DYLD_LIBRARY_PATH}"
export ROOT_INCLUDE_PATH="$(brew --prefix boost)/include:${ROOT_INCLUDE_PATH}"

# Run the executable with all arguments passed to this script
"${SCRIPT_DIR}/bin/iterateKT" "$@"
EOF

# Make run script executable
chmod +x ${ITERATEKT}/run_iterateKT.sh

echo "==============================================================="
echo "Build and setup completed successfully!"
echo "To run iterateKT, use:"
echo "  ./run_iterateKT.sh scripts/eta_decay.cpp"
echo ""
echo "Or for any other example:"
echo "  ./run_iterateKT.sh <path/to/script.cpp>"
echo "==============================================================="