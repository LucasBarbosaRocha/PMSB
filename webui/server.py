#!/usr/bin/env python3
"""
Servidor local (só biblioteca padrão do Python) para compilar e executar os
binários do PMSB pela web: escolher a ferramenta, subir grafo.fasta e
sequence.fasta, apertar "Executar" e ver o resultado na tela. Cada execução
é salva em webui/results/.

Uso:
    python3 webui/server.py
    (abre em http://localhost:8765)

Não precisa instalar nada — usa só a biblioteca padrão do Python 3.
"""
import cgi
import io
import json
import os
import re
import subprocess
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEBUI_DIR = os.path.dirname(os.path.abspath(__file__))
PORT = 8765

# ---------------------------------------------------------------------------
# Catálogo de ferramentas. Cada entrada aponta pra uma pasta do repositório
# (relativa à raiz do PMSB), o(s) comando(s) de compilação, o binário
# resultante e quais parâmetros de linha de comando ela espera.
#
# bmt_h1_b usa a biblioteca Bifrost (utils/myBifrost/bifrost/), que precisa
# estar compilada antes (cmake + make em utils/myBifrost/bifrost/build/) —
# se não estiver, o botão Compilar avisa e não quebra o resto do site. O
# compile usa -march=native de propósito: a lib Bifrost é compilada com essa
# flag (ativa AVX2), e alguns headers dela têm campos de struct condicionados
# a "#if defined(__AVX2__)" — compilar sem a mesma flag gera um binário que
# linka mas CRASHA em runtime por corrupção de memória (ODR violation), não
# por erro de compilação. Fica registrado aqui pra não se perder de novo.
#
# Fica de fora (dependência externa não resolvida por um "compilar com um
# clique"): variante2/btn_b — a pasta variante2/utils/ nem tem o wrapper
# myBifrost/ copiado. O protótipo variante2/btn (sem custo/emparelhamento,
# não implementa o algoritmo da tese) foi removido do catálogo por decisão
# do usuário — só btl (BMTC) representa a variante 2 aqui.
# ---------------------------------------------------------------------------
TOOLS = {
    "bmt": {
        "label": "Variante 1 (algoritmo exato) — BSMT: custo e mapeamento",
        "dir": "variante1/bmt",
        "binary": "bin/bmt",
        "compile": [["g++", "-O2", "-std=c++17", "-o", "bin/bmt", "bmt.cpp"]],
        "args": "full",  # -s -g -k -t -p
        "cost_regex": r"Custo:\s*(-?\d+)",
    },
    "oc_bmt": {
        "label": "Variante 1 (algoritmo exato) — BSMTd: apenas custo",
        "dir": "variante1/bmt_only_cost",
        "binary": "bin/bmt",
        "compile": [
            [
                "g++", "-O2", "-march=native", "-std=c++17", "-w",
                "-o", "bin/bmt", "bmt.cpp",
                "-I", "../utils/myBifrost/bifrost/src",
                "../utils/myBifrost/bifrost/build/src/libbifrost.a",
                "-lpthread", "-lz",
            ]
        ],
        "requires": "../utils/myBifrost/bifrost/build/src/libbifrost.a",
        "args": "full",
        "cost_regex": r"Custo:\s*(-?\d+)",
    },
    "bmt_h1": {
        "label": "Variante 1 (Heurística 1) — BSMTh1: custo e mapeamento",
        "dir": "variante1/bmt",
        "binary": "bin/bmt_h1",
        "compile": [["g++", "-O2", "-std=c++17", "-o", "bin/bmt_h1", "bmt_h1.cpp"]],
        "args": "full",
        "cost_regex": r"Custo:\s*(-?\d+)",
    },
    "bmt_h1_b": {
        "label": "Variante 1 (Heurística 1 com Bifrost) — BSMTh1b: custo e mapeamento",
        "dir": "variante1/bmt",
        "binary": "bin/bmt_h1_b",
        "compile": [
            [
                "g++", "-O2", "-march=native", "-std=c++17", "-w",
                "-o", "bin/bmt_h1_b", "bmt_h1_b.cpp",
                "-I", "../utils/myBifrost/bifrost/src",
                "../utils/myBifrost/bifrost/build/src/libbifrost.a",
                "-lpthread", "-lz",
            ]
        ],
        "requires": "../utils/myBifrost/bifrost/build/src/libbifrost.a",
        "args": "full",
        "cost_regex": r"Custo:\s*(-?\d+)",
    },
    "bmt_h2": {
        "label": "Variante 1 (Heurística 2) — BSMTh2: custo e mapeamento",
        "dir": "variante1/bmt",
        "binary": "bin/bmt_h2",
        "compile": [["g++", "-O2", "-std=c++17", "-o", "bin/bmt_h2", "bmt_h2.cpp"]],
        "args": "full",
        "cost_regex": r"Custo:\s*(-?\d+)",
    },
    "bmt_h3": {
        "label": "Variante 1 (Heurística 3) — BSMTh3: custo e mapeamento",
        "dir": "variante1/bmt",
        "binary": "bin/bmt_h3",
        "compile": [["g++", "-O2", "-std=c++17", "-o", "bin/bmt_h3", "bmt_h3.cpp"]],
        "args": "full",
        "cost_regex": r"Custo:\s*(-?\d+)",
    },
    "oc_bmt_h1": {
        "label": "Variante 1 (Heurística 1) — BSMTh1d: apenas custo",
        "dir": "variante1/bmt_only_cost",
        "binary": "bin/bmt_h1",
        "compile": [
            [
                "g++", "-O2", "-march=native", "-std=c++17", "-w",
                "-o", "bin/bmt_h1", "bmt_h1.cpp",
                "-I", "../utils/myBifrost/bifrost/src",
                "../utils/myBifrost/bifrost/build/src/libbifrost.a",
                "-lpthread", "-lz",
            ]
        ],
        "requires": "../utils/myBifrost/bifrost/build/src/libbifrost.a",
        "args": "full",
        "cost_regex": r"Custo:\s*(-?\d+)",
    },
    "oc_bmt_h1_b": {
        "label": "Variante 1 (Heurística 1 com Bifrost) — BSMTh1db: apenas custo",
        "dir": "variante1/bmt_only_cost",
        "binary": "bin/bmt_h1_b",
        "compile": [
            [
                "g++", "-O2", "-march=native", "-std=c++17", "-w",
                "-o", "bin/bmt_h1_b", "bmt_h1_b.cpp",
                "-I", "../utils/myBifrost/bifrost/src",
                "../utils/myBifrost/bifrost/build/src/libbifrost.a",
                "-lpthread", "-lz",
            ]
        ],
        "requires": "../utils/myBifrost/bifrost/build/src/libbifrost.a",
        "args": "full",
        "cost_regex": r"Custo:\s*(-?\d+)",
    },
    "btl": {
        "label": "Variante 2 — BMTC (emparelhamento Húngaro)",
        "dir": "variante2/btl",
        "binary": "btl",
        "compile": [
            ["g++", "-O2", "-std=c++17", "-c", "btl.cpp", "-o", "main.o"],
            ["g++", "-O2", "-std=c++17", "-c", "Hungarian.cpp", "-o", "hung.o"],
            ["g++", "-o", "btl", "main.o", "hung.o"],
        ],
        "args": "basic",
        "cost_regex": r"Custo:\s*(-?\d+(?:\.\d+)?)",
    },
    "pesb": {
        "label": "Auxiliar — Pesb: k-mers da sequência ausentes no grafo",
        "dir": "variante1/bmt",
        "binary": "bin/pesb",
        "compile": [["g++", "-O2", "-std=c++17", "-o", "bin/pesb", "pesb.cpp"]],
        "args": "full",
        "cost_regex": r"Vamos inserir:\s*(-?\d+)",
    },
}


def tool_dir(tool):
    return os.path.join(ROOT, tool["dir"])


def safe_filename(name):
    name = os.path.basename(name or "arquivo")
    name = re.sub(r"[^A-Za-z0-9_.\-]", "_", name)
    return name or "arquivo"


def run_compile(tool_id):
    tool = TOOLS[tool_id]
    cwd = tool_dir(tool)
    os.makedirs(os.path.join(cwd, "bin"), exist_ok=True)
    os.makedirs(os.path.join(cwd, "uploads"), exist_ok=True)

    requires = tool.get("requires")
    if requires and not os.path.isfile(os.path.join(cwd, requires)):
        return False, (
            f"Dependência ausente: {requires}\n"
            "A biblioteca Bifrost ainda não foi compilada. Rode:\n"
            "  cd variante1/utils/myBifrost/bifrost\n"
            "  mkdir -p build && cd build && cmake .. && make -j$(nproc)"
        )

    log = []
    ok = True
    for cmd in tool["compile"]:
        log.append("$ " + " ".join(cmd))
        try:
            proc = subprocess.run(
                cmd, cwd=cwd, capture_output=True, text=True, timeout=180
            )
        except subprocess.TimeoutExpired:
            log.append("[erro] compilação excedeu o tempo limite (180s)")
            ok = False
            break
        if proc.stdout:
            log.append(proc.stdout.rstrip())
        if proc.stderr:
            log.append(proc.stderr.rstrip())
        if proc.returncode != 0:
            ok = False
            break
    return ok, "\n".join(log)


def run_tool(tool_id, k, graph_bytes, graph_name, seq_bytes, seq_name):
    # Grafo sempre tradicional (-t 0) e mapeamento sempre por passeio (-p 0) —
    # não é exposto na interface, fixo por decisão do usuário.
    t = 0
    p = 0
    tool = TOOLS[tool_id]
    cwd = tool_dir(tool)
    binary = os.path.join(cwd, tool["binary"])
    if not os.path.isfile(binary):
        return {
            "ok": False,
            "error": f"Binário não encontrado ({tool['binary']}). Compile primeiro.",
        }

    stamp = time.strftime("%Y%m%d_%H%M%S")
    uploads_dir = os.path.join(cwd, "uploads")
    os.makedirs(uploads_dir, exist_ok=True)

    graph_path = os.path.join(uploads_dir, f"{stamp}_{safe_filename(graph_name)}")
    seq_path = os.path.join(uploads_dir, f"{stamp}_{safe_filename(seq_name)}")
    with open(graph_path, "wb") as f:
        f.write(graph_bytes)
    with open(seq_path, "wb") as f:
        f.write(seq_bytes)

    args = [binary, "-s", seq_path, "-g", graph_path, "-k", str(k), "-t", str(t)]
    if tool["args"] == "full":
        args += ["-p", str(p)]

    try:
        proc = subprocess.run(args, cwd=cwd, capture_output=True, text=True, timeout=60)
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": "Execução excedeu o tempo limite (60s)."}

    stdout = proc.stdout
    stderr = proc.stderr
    cost = None
    if tool["cost_regex"]:
        matches = re.findall(tool["cost_regex"], stdout)
        if matches:
            cost = matches[-1]  # última ocorrência (útil quando a saída não tem rótulo "Cost:")

    results_dir = os.path.join(WEBUI_DIR, "results")
    os.makedirs(results_dir, exist_ok=True)
    result_name = f"{stamp}_{tool_id}.txt"
    result_path = os.path.join(results_dir, result_name)
    with open(result_path, "w") as f:
        f.write(f"ferramenta: {tool_id}\n")
        f.write(f"data: {stamp}\n")
        f.write(f"grafo: {graph_name}\n")
        f.write(f"sequencia: {seq_name}\n")
        f.write(f"k={k} t={t}" + (f" p={p}" if tool['args'] == 'full' else "") + "\n")
        f.write(f"comando: {' '.join(args)}\n")
        f.write("\n--- stdout ---\n")
        f.write(stdout)
        if stderr:
            f.write("\n--- stderr ---\n")
            f.write(stderr)

    return {
        "ok": True,
        "returncode": proc.returncode,
        "stdout": stdout,
        "stderr": stderr,
        "cost": cost,
        "saved_to": os.path.relpath(result_path, ROOT),
    }


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        pass  # silencia log padrão barulhento

    def _send_json(self, obj, status=200):
        body = json.dumps(obj).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_file(self, path, content_type):
        try:
            with open(path, "rb") as f:
                body = f.read()
        except FileNotFoundError:
            self.send_response(404)
            self.end_headers()
            return
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path == "/" or self.path == "/index.html":
            self._send_file(os.path.join(WEBUI_DIR, "index.html"), "text/html; charset=utf-8")
        elif self.path == "/app.js":
            self._send_file(os.path.join(WEBUI_DIR, "app.js"), "application/javascript; charset=utf-8")
        elif self.path == "/style.css":
            self._send_file(os.path.join(WEBUI_DIR, "style.css"), "text/css; charset=utf-8")
        elif self.path == "/api/tools":
            payload = [
                {"id": tid, "label": t["label"], "args": t["args"], "hasCost": bool(t["cost_regex"])}
                for tid, t in TOOLS.items()
            ]
            self._send_json(payload)
        else:
            self.send_response(404)
            self.end_headers()

    def do_POST(self):
        if self.path == "/api/compile":
            length = int(self.headers.get("Content-Length", 0))
            data = json.loads(self.rfile.read(length) or b"{}")
            tool_id = data.get("tool")
            if tool_id not in TOOLS:
                return self._send_json({"ok": False, "error": "ferramenta desconhecida"}, 400)
            ok, log = run_compile(tool_id)
            self._send_json({"ok": ok, "log": log})
            return

        if self.path == "/api/run":
            ctype, pdict = cgi.parse_header(self.headers.get("Content-Type", ""))
            if ctype != "multipart/form-data":
                return self._send_json({"ok": False, "error": "esperado multipart/form-data"}, 400)
            pdict["boundary"] = pdict["boundary"].encode("utf-8")
            length = int(self.headers.get("Content-Length", 0))
            form = cgi.parse_multipart(io.BytesIO(self.rfile.read(length)), pdict)

            def field(name, default=""):
                v = form.get(name)
                if not v:
                    return default
                v0 = v[0]
                return v0.decode("utf-8") if isinstance(v0, bytes) else v0

            tool_id = field("tool")
            if tool_id not in TOOLS:
                return self._send_json({"ok": False, "error": "ferramenta desconhecida"}, 400)

            k = field("k", "3")

            graph_list = form.get("graph")
            seq_list = form.get("sequence")
            if not graph_list or not seq_list:
                return self._send_json({"ok": False, "error": "envie os dois arquivos (grafo e sequência)"}, 400)

            graph_bytes = graph_list[0] if isinstance(graph_list[0], bytes) else graph_list[0].encode("utf-8")
            seq_bytes = seq_list[0] if isinstance(seq_list[0], bytes) else seq_list[0].encode("utf-8")

            graph_name = field("graph_name", "grafo.fasta")
            seq_name = field("sequence_name", "sequence.fasta")

            try:
                k_int = int(k)
            except ValueError:
                return self._send_json({"ok": False, "error": "k precisa ser inteiro"}, 400)

            result = run_tool(tool_id, k_int, graph_bytes, graph_name, seq_bytes, seq_name)
            self._send_json(result)
            return

        self.send_response(404)
        self.end_headers()


def main():
    server = ThreadingHTTPServer(("localhost", PORT), Handler)
    print(f"PMSB webui rodando em http://localhost:{PORT}  (Ctrl+C pra parar)")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
