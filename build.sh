#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
if ! command -v c++ >/dev/null 2>&1; then
  echo '需要 C++ 编译工具 / C++ compiler required. On Mac run: xcode-select --install'
  exit 1
fi
mkdir -p build
c++ -std=c++17 -O2 -Wall -Wextra -Wpedantic src/main.cpp -o build/peptide-lens
printf '编译完成 / Built: build/peptide-lens\n'
