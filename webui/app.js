let TOOLS = [];

async function loadTools() {
  const res = await fetch("/api/tools");
  TOOLS = await res.json();
  const sel = document.getElementById("tool");
  sel.innerHTML = "";
  for (const t of TOOLS) {
    const opt = document.createElement("option");
    opt.value = t.id;
    opt.textContent = t.label;
    sel.appendChild(opt);
  }
  updateHint();
}

function currentTool() {
  const id = document.getElementById("tool").value;
  return TOOLS.find((t) => t.id === id);
}

function updateHint() {
  const tool = currentTool();
  if (!tool) return;
  document.getElementById("toolHint").textContent = tool.hasCost
    ? "Esta ferramenta reporta um custo — aparecerá em destaque no resultado."
    : "Esta ferramenta não reporta um custo único (veja a saída completa).";
}

document.addEventListener("DOMContentLoaded", () => {
  loadTools();
  document.getElementById("tool").addEventListener("change", updateHint);

  document.getElementById("btnCompile").addEventListener("click", async () => {
    const tool = currentTool();
    const log = document.getElementById("compileLog");
    log.hidden = false;
    log.textContent = "Compilando " + tool.id + " ...";
    const res = await fetch("/api/compile", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ tool: tool.id }),
    });
    const data = await res.json();
    log.textContent = (data.ok ? "OK\n\n" : "FALHOU\n\n") + (data.log || data.error || "");
  });

  document.getElementById("btnRun").addEventListener("click", async () => {
    const tool = currentTool();
    const graphInput = document.getElementById("graph");
    const seqInput = document.getElementById("sequence");
    const output = document.getElementById("output");
    const costLine = document.getElementById("costLine");
    const savedLine = document.getElementById("savedLine");

    if (!graphInput.files[0] || !seqInput.files[0]) {
      output.textContent = "Envie o arquivo do grafo e o da sequência antes de executar.";
      costLine.hidden = true;
      return;
    }

    const form = new FormData();
    form.append("tool", tool.id);
    form.append("k", document.getElementById("k").value);
    form.append("graph", graphInput.files[0]);
    form.append("graph_name", graphInput.files[0].name);
    form.append("sequence", seqInput.files[0]);
    form.append("sequence_name", seqInput.files[0].name);

    output.textContent = "Executando ...";
    costLine.hidden = true;
    savedLine.textContent = "";

    const res = await fetch("/api/run", { method: "POST", body: form });
    const data = await res.json();

    if (!data.ok) {
      output.textContent = "Erro: " + (data.error || "desconhecido");
      return;
    }

    output.textContent = data.stdout + (data.stderr ? "\n--- stderr ---\n" + data.stderr : "");
    if (data.cost !== null && data.cost !== undefined) {
      costLine.hidden = false;
      costLine.textContent = "Custo: " + data.cost;
    } else {
      costLine.hidden = true;
    }
    savedLine.textContent = "Resultado salvo em: " + data.saved_to;
  });
});
