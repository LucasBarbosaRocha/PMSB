#!/usr/bin/env bash
# Compila todos os programas da pasta bmt/ e coloca os binários em bmt/bin/.
#
# Uso:
#   ./build.sh
#
# Os binários gerados usam includes relativos ("../utils/...cpp"), então este
# script precisa ser executado a partir da pasta bmt/ (é o que ele já faz,
# usando o diretório onde o próprio script está).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

BIN_DIR="$SCRIPT_DIR/bin"
mkdir -p "$BIN_DIR"

CXX=${CXX:-g++}
CXXFLAGS="-O2 -std=c++17"

# Programas que dependem apenas de utils/*.cpp (sem dependências externas).
PROGRAMS=(bmt bmt_h1 bmt_h2 bmt_h3 pesb)

echo "Compilando com: $CXX $CXXFLAGS"
echo

for prog in "${PROGRAMS[@]}"; do
    echo "-> $prog"
    "$CXX" $CXXFLAGS -o "$BIN_DIR/$prog" "$prog.cpp"
done

echo
echo "Programas compilados em: $BIN_DIR"
echo

# bmt_h1_b depende da biblioteca Bifrost (utils/myBifrost/bifrost), que não
# vem pré-compilada neste repositório. Compile-a primeiro:
#   cd utils/myBifrost/bifrost && mkdir -p build && cd build && cmake .. && make -j
#
# IMPORTANTE sobre -march=native: o CMakeLists da Bifrost compila a lib com
# COMPILATION_ARCH=native (ativa AVX2 nesta máquina). Vários headers da
# Bifrost (ex. BlockedBloomFilter.hpp) têm campos de struct condicionados a
# "#if defined(__AVX2__)" — se bmt_h1_b.cpp for compilado SEM -march=native/
# -mavx2, o tamanho dessas structs diverge do que está na lib já compilada
# (violação de ODR) e o binário CRASHA em runtime com corrupção de memória
# silenciosa (não é erro de link nem de compilação — só aparece rodando).
# Por isso usamos -march=native aqui também, sempre a mesma flag da lib.
BIFROST_DIR="../utils/myBifrost/bifrost"
BIFROST_LIB="$BIFROST_DIR/build/src/libbifrost.a"

echo "-> bmt_h1_b (requer Bifrost já compilada em $BIFROST_DIR/build)"
if [ -f "$BIFROST_LIB" ]; then
    if "$CXX" -O2 -march=native -std=c++17 -o "$BIN_DIR/bmt_h1_b" bmt_h1_b.cpp \
        -I "$BIFROST_DIR/src" "$BIFROST_LIB" -lpthread -lz 2>"$BIN_DIR/bmt_h1_b_build.log"; then
        echo "   OK: $BIN_DIR/bmt_h1_b"
    else
        echo "   AVISO: falha ao compilar/linkar bmt_h1_b."
        echo "   Detalhes em: $BIN_DIR/bmt_h1_b_build.log"
    fi
else
    echo "   AVISO: Bifrost ainda não compilada ($BIFROST_LIB não existe)."
    echo "   Rode: cd $BIFROST_DIR && mkdir -p build && cd build && cmake .. && make -j\$(nproc)"
fi

echo
echo "Build concluído."
