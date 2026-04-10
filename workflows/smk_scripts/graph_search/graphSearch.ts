export type VisNode = {
  id?: string | number;
  label?: string;
  title?: string;
};

export type VisNodesDataSet = {
  get(): VisNode[];
};

export type VisNetwork = {
  getConnectedNodes(nodeId: string | number): Array<string | number>;
  selectNodes(nodeIds: Array<string | number>, highlightEdges?: boolean): void;
  focus(nodeId: string | number, options?: unknown): void;
  fit(options?: unknown): void;
  unselectAll(): void;
};

export function getNodeText(n: VisNode): string {
  let t = "";
  if (n.id !== undefined && n.id !== null) t += String(n.id) + "\n";
  if (n.label) t += String(n.label) + "\n";
  if (n.title) t += String(n.title) + "\n";
  return t;
}

export function compileRegex(pattern: string, flags = "i"): RegExp {
  const raw = (pattern || "").trim();
  if (!raw) throw new Error("empty regex");
  return new RegExp(raw, flags);
}

export function findMatchingNodeIds(nodes: VisNode[], re: RegExp): Array<string | number> {
  const out: Array<string | number> = [];
  for (const n of nodes) {
    const id = n.id;
    if (id === undefined || id === null) continue;
    if (re.test(getNodeText(n))) out.push(id);
  }
  return out;
}

export function expandWithNeighbors(network: VisNetwork, matched: Array<string | number>): Array<string | number> {
  const keep = new Set<string | number>();
  for (const id of matched) {
    keep.add(id);
    const neigh = network.getConnectedNodes(id) || [];
    for (const nid of neigh) keep.add(nid);
  }
  return Array.from(keep);
}

export async function waitFor<T>(
  resolve: () => T | null | undefined,
  opts: { timeoutMs: number; intervalMs: number }
): Promise<T> {
  const start = Date.now();
  // eslint-disable-next-line no-constant-condition
  while (true) {
    const v = resolve();
    if (v) return v;
    if (Date.now() - start >= opts.timeoutMs) throw new Error("timeout");
    await new Promise((r) => setTimeout(r, opts.intervalMs));
  }
}

export type InitGraphSearchOptions = {
  container: HTMLElement;
  resolveNetwork: () => VisNetwork | null | undefined;
  resolveNodesDS: (network: VisNetwork) => VisNodesDataSet | null | undefined;
  placeholder?: string;
  onStatus?: (msg: string) => void;
  readyTimeoutMs?: number;
  readyIntervalMs?: number;
};

export async function initGraphSearch(opts: InitGraphSearchOptions): Promise<void> {
  const readyTimeoutMs = opts.readyTimeoutMs ?? 8000;
  const readyIntervalMs = opts.readyIntervalMs ?? 75;

  const panel = document.createElement("div");
  panel.style.position = "absolute";
  panel.style.top = "10px";
  panel.style.left = "10px";
  panel.style.zIndex = "20";
  panel.style.background = "rgba(255,255,255,0.92)";
  panel.style.padding = "8px 10px";
  panel.style.border = "1px solid #ddd";
  panel.style.borderRadius = "6px";
  panel.style.fontFamily = "sans-serif";
  panel.style.fontSize = "12px";
  panel.style.maxWidth = "520px";
  panel.innerHTML = `
    <div style="margin-bottom:6px;">
      <b>Search</b> (JS regex; matches node id/tooltip). Selects matches + 1-hop neighbors.
    </div>
    <div style="display:flex; gap:6px; align-items:center;">
      <input data-reportr="regex" type="text" placeholder="${opts.placeholder ?? ""}" style="flex:1; padding:4px 6px;"/>
      <button data-reportr="apply" style="padding:4px 8px;" disabled>Apply</button>
      <button data-reportr="clear" style="padding:4px 8px;" disabled>Clear</button>
    </div>
    <div data-reportr="msg" style="margin-top:6px; color:#444;"></div>
  `;

  opts.container.style.position = "relative";
  opts.container.appendChild(panel);

  const msgEl = panel.querySelector('[data-reportr="msg"]') as HTMLDivElement;
  const inputEl = panel.querySelector('[data-reportr="regex"]') as HTMLInputElement;
  const applyBtn = panel.querySelector('[data-reportr="apply"]') as HTMLButtonElement;
  const clearBtn = panel.querySelector('[data-reportr="clear"]') as HTMLButtonElement;

  function setMsg(s: string) {
    msgEl.textContent = s;
    opts.onStatus?.(s);
  }

  let network: VisNetwork | null = null;
  let nodesDS: VisNodesDataSet | null = null;

  setMsg("Initializing graph…");
  try {
    network = await waitFor(() => opts.resolveNetwork(), {
      timeoutMs: readyTimeoutMs,
      intervalMs: readyIntervalMs,
    });
    nodesDS = await waitFor(() => opts.resolveNodesDS(network as VisNetwork), {
      timeoutMs: readyTimeoutMs,
      intervalMs: readyIntervalMs,
    });
  } catch {
    setMsg("Graph not ready yet (network not initialized). Reload the page or try again.");
    return;
  }

  applyBtn.disabled = false;
  clearBtn.disabled = false;
  setMsg("");

  function applyRegex() {
    const raw = (inputEl.value || "").trim();
    if (!raw) {
      setMsg("Enter a regex.");
      return;
    }

    let re: RegExp;
    try {
      re = compileRegex(raw, "i");
    } catch (e) {
      setMsg("Invalid regex: " + String(e));
      return;
    }

    const nodes = (nodesDS as VisNodesDataSet).get();
    const matched = findMatchingNodeIds(nodes, re);
    if (matched.length === 0) {
      setMsg("No matches.");
      return;
    }

    const ids = expandWithNeighbors(network as VisNetwork, matched);
    (network as VisNetwork).selectNodes(ids, false);
    if (matched.length === 1) {
      (network as VisNetwork).focus(matched[0], { scale: 1.2, animation: { duration: 300 } });
    } else {
      (network as VisNetwork).fit({ nodes: ids, animation: { duration: 300 } });
    }
    setMsg(`Matched ${matched.length} node(s); selected ${ids.length} incl. neighbors.`);
  }

  applyBtn.addEventListener("click", applyRegex);
  inputEl.addEventListener("keydown", (e) => {
    if ((e as KeyboardEvent).key === "Enter") applyRegex();
  });
  clearBtn.addEventListener("click", () => {
    inputEl.value = "";
    setMsg("");
    (network as VisNetwork).unselectAll();
    (network as VisNetwork).fit({ animation: { duration: 200 } });
  });
}

