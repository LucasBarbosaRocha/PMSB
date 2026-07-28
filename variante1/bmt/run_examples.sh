#!/usr/bin/env bash
# Roda os programas de bmt/ com os dados de exemplo em examples/.
# Compila antes, se os binários ainda não existirem.
#
# Uso:
#   ./run_examples.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

BIN_DIR="$SCRIPT_DIR/bin"
EXAMPLES_DIR="$SCRIPT_DIR/examples"

if [ ! -x "$BIN_DIR/bmt" ]; then
    echo "Binários não encontrados, compilando primeiro..."
    ./build.sh
    echo
fi

GRAFO="$EXAMPLES_DIR/grafo.fasta"
SEQ="$EXAMPLES_DIR/sequence.fasta"
K=3

run() {
    local titulo="$1"; shift
    echo "=============================================="
    echo "$titulo"
    echo "\$ $*"
    echo "----------------------------------------------"
    "$@"
    echo
}

run "bmt (grafo tradicional, -t 0)" \
    "$BIN_DIR/bmt" -s "$SEQ" -g "$GRAFO" -k "$K" -t 0 -p 0

run "bmt (grafo simplificado, -t 1)" \
    "$BIN_DIR/bmt" -s "$SEQ" -g "$GRAFO" -k "$K" -t 1 -p 0

run "bmt_h1 (heurística 1, mapeamento por âncoras)" \
    "$BIN_DIR/bmt_h1" -s "$SEQ" -g "$GRAFO" -k "$K" -t 0 -p 0

run "bmt_h2 (heurística seed-and-extend)" \
    "$BIN_DIR/bmt_h2" -s "$SEQ" -g "$GRAFO" -k "$K" -t 0 -p 0

run "bmt_h3 (seed-and-extend com posições válidas)" \
    "$BIN_DIR/bmt_h3" -s "$SEQ" -g "$GRAFO" -k "$K" -t 0 -p 0

run "pesb (contagem de k-mers ausentes no grafo)" \
    "$BIN_DIR/pesb" -s "$SEQ" -g "$GRAFO" -k "$K" -t 0 -p 0

echo "=============================================="
echo "Todos os exemplos executados."
