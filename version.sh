#!/usr/bin/env bash

set -euo pipefail

# ---- Configuration ----
CMAKE_FILE="CMakeLists.txt"

# ---- Input validation ----
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 X.Y.Z"
    exit 1
fi

VERSION="$1"

# Validate semantic version X.Y.Z where X,Y,Z are non-negative integers
if ! [[ "$VERSION" =~ ^([0-9]+)\.([0-9]+)\.([0-9]+)$ ]]; then
    echo "Error: Version must be in format X.Y.Z (numeric)."
    exit 1
fi

# ---- Sanity checks ----
if [ ! -f "$CMAKE_FILE" ]; then
    echo "Error: CMake file not found: $CMAKE_FILE"
    exit 1
fi

# ---- Update VERSIONs in C++ file ----
# This replaces the numeric value after the define
sed -i.bak -E "s/(PROJECT\(TMDET VERSION )[0-9]+\.[0-9]+\.[0-9]+\)/\1$VERSION)/" "$CMAKE_FILE"

rm -f "${CMAKE_FILE}.bak"

echo "Updated to $VERSION in $CMAKE_FILE"

# ---- Git checks ----
if ! git diff --quiet; then
    git add "$CMAKE_FILE"
    git commit -m "chore: bump version to $VERSION"
else
    echo "No changes detected in version file."
fi

TAG="v$VERSION"

# Ensure tag does not already exist
if git rev-parse "$TAG" >/dev/null 2>&1; then
    echo "Error: Git tag $TAG already exists."
    exit 1
fi

# ---- Create git tag ----
git tag -a "$TAG"
git push origin "$TAG"
echo "Created git tag: $TAG"
echo "Version synchronization complete."
