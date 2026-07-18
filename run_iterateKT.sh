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
