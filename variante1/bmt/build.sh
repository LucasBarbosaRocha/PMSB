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
# vem pré-compilada neste repositório. A compilação (checagem de sintaxe)
# funciona, mas o link só funciona se a libbifrost já tiver sido construída
# (veja utils/myBifrost/bifrost/CMakeLists.txt). Tentamos compilar e linkar;
# se falhar, avisamos e seguimos sem interromper o build dos demais.
echo "-> bmt_h1_b (requer Bifrost já compilado; pode falhar)"
if "$CXX" $CXXFLAGS -o "$BIN_DIR/bmt_h1_b" bmt_h1_b.cpp -lz -lpthread 2>"$BIN_DIR/bmt_h1_b_build.log"; then
    echo "   OK: $BIN_DIR/bmt_h1_b"
else
    echo "   AVISO: não foi possível linkar bmt_h1_b (biblioteca Bifrost ausente)."
    echo "   Detalhes em: $BIN_DIR/bmt_h1_b_build.log"
    echo "   Para habilitar, compile a Bifrost em utils/myBifrost/bifrost/ antes de rodar este script."
fi

echo
echo "Build concluído."
