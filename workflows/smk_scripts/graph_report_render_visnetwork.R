#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(visNetwork)
  library(htmlwidgets)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: graph_report_render_visnetwork.R <nodes.tsv> <edges.full.tsv> <clusters.tsv> <cluster_annotations.tsv> <out.html> <out.nodes.enriched.tsv>")
}

nodes_path <- args[[1]]
edges_path <- args[[2]]
clusters_path <- args[[3]]
ann_path <- args[[4]]
out_html <- args[[5]]
out_nodes_enriched <- args[[6]]

nodes <- read_tsv(nodes_path, show_col_types = FALSE)
edges <- read_tsv(edges_path, show_col_types = FALSE)
clusters <- read_tsv(clusters_path, show_col_types = FALSE)
ann <- read_tsv(ann_path, show_col_types = FALSE)

nodes2 <- nodes %>%
  left_join(clusters, by = c("repeat" = "repeat")) %>%
  mutate(cluster_id = ifelse(is.na(cluster_id), "UNCLUSTERED", cluster_id))

ann2 <- ann %>%
  select(any_of(c(
    "cluster_id", "n_repeats", "best_hit", "best_species", "best_title",
    "best_name", "best_pident", "best_evalue", "oligo_best_sseqid", "oligo_fitting"
  ))) %>%
  mutate(across(everything(), ~ifelse(is.na(.x), "", as.character(.x))))

norm_id <- function(x) {
  # LocalDB subjects are stored as "<project>|<subject_type>|<subject_id>|<repeat_id>".
  # The graph nodes use "<repeat_id>" (e.g. "SMPL=...|..."). Normalize BLAST ids.
  x2 <- as.character(x)
  x2 <- ifelse(str_detect(x2, "\\|SMPL="), sub(".*\\|(SMPL=.*)$", "\\1", x2), x2)
  x2
}

nodes3 <- nodes2 %>%
  left_join(ann2, by = "cluster_id") %>%
  mutate(
    id = .data[["repeat"]],
    # Labels are hidden by default to avoid overflow. The UI can turn them on.
    label = "",
    label_repeat = .data[["repeat"]],
    title = "",
    group = .data[["cluster_id"]],
    shape = ifelse(.data[["subject_type"]] == "comparative", "triangle", "dot"),
    # Use satMiner abundance (% of repeatome) to scale node size.
    value = ifelse(is.na(.data[["size_pct"]]) | .data[["size_pct"]] <= 0, 1, .data[["size_pct"]])
  )

edges2 <- edges %>%
  transmute(
    from = norm_id(qseqid),
    to = norm_id(sseqid),
    value = pmax(coverage, 0),
    coverage = coverage,
    coverage_type = coverage_type,
    pident = pident,
    evalue = evalue,
    task = task,
    title = paste0(
      "coverage_type=", coverage_type,
      " cov=", round(coverage, 3),
      " pident=", round(pident, 1),
      " task=", task
    )
  ) %>%
  filter(from %in% nodes3$id, to %in% nodes3$id)

dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
write_tsv(nodes3, out_nodes_enriched)

g <- visNetwork(nodes3, edges2, main = "") %>%
  visNodes(font = list(size = 10), scaling = list(label = list(enabled = TRUE, min = 8, max = 18))) %>%
  visEdges(smooth = FALSE, color = list(opacity = 0.35)) %>%
  visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = TRUE) %>%
  visLayout(randomSeed = 1) %>%
  onRender("
function(el, x) {
  var container = document.getElementById(el.id);
  if (!container) return;
  container.style.position = 'relative';
  // Make the widget fill the viewport (avoid fixed 960x500 box sizing).
  try {
    document.documentElement.style.height = '100%';
    document.body.style.height = '100%';
    document.body.style.margin = '0';
  } catch (e) {}
  try {
    el.style.width = '100vw';
    el.style.height = '100vh';
  } catch (e) {}
  try {
    var graphEl0 = document.getElementById('graph' + el.id);
    if (graphEl0) { graphEl0.style.width = '100%'; graphEl0.style.height = '100%'; }
  } catch (e) {}

  function _waitFor(resolve, timeoutMs, intervalMs) {
    var start = Date.now();
    return new Promise(function(ok, bad) {
      function tick() {
        try {
          var v = resolve();
          if (v) return ok(v);
        } catch (e) {}
        if (Date.now() - start >= timeoutMs) return bad(new Error('timeout'));
        setTimeout(tick, intervalMs);
      }
      tick();
    });
  }

  function _getNodeText(n) {
    var t = '';
    if (n && n.id !== undefined && n.id !== null) t += String(n.id) + '\\n';
    if (n && n.label) t += String(n.label) + '\\n';
    if (n && n.title) t += String(n.title) + '\\n';
    return t;
  }

  function _compileRegex(raw) {
    raw = (raw || '').trim();
    if (!raw) throw new Error('empty regex');
    return new RegExp(raw, 'i');
  }

  function _findMatchingNodeIds(nodes, re) {
    var out = [];
    for (var i = 0; i < nodes.length; i++) {
      var n = nodes[i];
      if (!n || n.id === undefined || n.id === null) continue;
      if (re.test(_getNodeText(n))) out.push(n.id);
    }
    return out;
  }

  function _expandWithNeighbors(network, matched) {
    var keep = {};
    for (var i = 0; i < matched.length; i++) {
      var id = matched[i];
      keep[id] = true;
      var neigh = network.getConnectedNodes(id) || [];
      for (var j = 0; j < neigh.length; j++) keep[neigh[j]] = true;
    }
    return Object.keys(keep);
  }

  var panel = document.createElement('div');
  panel.style.position = 'absolute';
  panel.style.top = '10px';
  panel.style.right = '10px';
  panel.style.zIndex = 20;
  panel.style.background = 'rgba(255,255,255,0.92)';
  panel.style.padding = '8px 10px';
  panel.style.border = '1px solid #ddd';
  panel.style.borderRadius = '6px';
  panel.style.fontFamily = 'sans-serif';
  panel.style.fontSize = '12px';
  panel.style.maxWidth = '520px';
  panel.innerHTML = '' +
    '<div style=\"margin-bottom:6px;\">' +
      '<b>Search</b> (JS regex; matches node id/tooltip). Selects matches + 1-hop neighbors.' +
    '</div>' +
    '<div style=\"display:flex; gap:6px; align-items:center;\">' +
      '<input data-reportr=\"regex\" type=\"text\" placeholder=\"e.g. ^SMPL=KA1\\\\b | ORG=Aegilops | pdiv_mean\" style=\"flex:1; padding:4px 6px;\"/>' +
      '<button data-reportr=\"apply\" style=\"padding:4px 8px;\" disabled>Apply</button>' +
      '<button data-reportr=\"clear\" style=\"padding:4px 8px;\" disabled>Clear</button>' +
    '</div>' +
    '<hr style=\"margin:8px 0; border:0; border-top:1px solid #eee;\"/>' +
    '<div style=\"display:grid; grid-template-columns: 1fr 1fr; gap:8px;\">' +
      '<div>' +
        '<div><b>Labels</b></div>' +
        '<label style=\"display:block; margin-top:4px;\"><input data-reportr=\"showLabels\" type=\"checkbox\"/> Show labels</label>' +
        '<label style=\"display:block; margin-top:4px;\">Mode ' +
          '<select data-reportr=\"labelMode\" style=\"width:100%; padding:2px 4px;\">' +
            '<option value=\"repeat\">repeat</option>' +
            '<option value=\"repeat_size\">repeat + size_pct</option>' +
            '<option value=\"repeat_pdiv\">repeat + pdiv_mean</option>' +
            '<option value=\"size_only\">size_pct only</option>' +
          '</select>' +
        '</label>' +
      '</div>' +
      '<div>' +
        '<div><b>Tooltip fields</b></div>' +
        '<label style=\"display:block; margin-top:4px;\"><input data-reportr=\"tt_size\" type=\"checkbox\" checked/> size_pct</label>' +
        '<label style=\"display:block; margin-top:4px;\"><input data-reportr=\"tt_pdiv\" type=\"checkbox\" checked/> pdiv_mean</label>' +
        '<label style=\"display:block; margin-top:4px;\"><input data-reportr=\"tt_tarean\" type=\"checkbox\" checked/> TAREAN_annotation</label>' +
      '</div>' +
    '</div>' +
    '<hr style=\"margin:8px 0; border:0; border-top:1px solid #eee;\"/>' +
    '<div style=\"display:grid; grid-template-columns: 1fr 1fr; gap:8px;\">' +
      '<div>' +
        '<div><b>Nodes</b></div>' +
        '<label style=\"display:block; margin-top:4px;\">Samples' +
          '<select data-reportr=\"samples\" multiple size=\"6\" style=\"width:100%; padding:2px 4px;\"></select>' +
        '</label>' +
        '<label style=\"display:block; margin-top:4px;\"><input data-reportr=\"keepSamples\" type=\"checkbox\" checked/> keep selected samples</label>' +
        '<label style=\"display:block; margin-top:4px;\"><input data-reportr=\"keepComparatives\" type=\"checkbox\" checked/> keep comparatives</label>' +
      '</div>' +
      '<div>' +
        '<div><b>Edges</b></div>' +
        '<label style=\"display:block; margin-top:4px;\">Tasks' +
          '<select data-reportr=\"tasks\" multiple size=\"6\" style=\"width:100%; padding:2px 4px;\"></select>' +
        '</label>' +
        '<label style=\"display:block; margin-top:4px;\">Min coverage_type ' +
          '<select data-reportr=\"minCovType\" style=\"width:100%; padding:2px 4px;\">' +
            '<option value=\"weak\">weak</option>' +
            '<option value=\"partial\">partial</option>' +
            '<option value=\"composite\">composite</option>' +
            '<option value=\"near-full\" selected>near-full</option>' +
          '</select>' +
        '</label>' +
        '<label style=\"display:block; margin-top:4px;\"><input data-reportr=\"hideIsolated\" type=\"checkbox\"/> hide isolated nodes</label>' +
      '</div>' +
    '</div>' +
    '<div style=\"display:flex; gap:6px; margin-top:8px;\">' +
      '<button data-reportr=\"applyFilters\" style=\"padding:4px 8px;\" disabled>Apply filters</button>' +
      '<button data-reportr=\"resetFilters\" style=\"padding:4px 8px;\" disabled>Reset filters</button>' +
    '</div>' +
    '<div data-reportr=\"msg\" style=\"margin-top:6px; color:#444;\"></div>';
  container.appendChild(panel);

  // Centered title overlay (visNetwork main text alignment is inconsistent across browsers/versions).
  var title = document.createElement('div');
  title.textContent = 'Repeat interaction graph';
  title.style.position = 'absolute';
  title.style.top = '10px';
  title.style.left = '50%';
  title.style.transform = 'translateX(-50%)';
  title.style.zIndex = 10;
  title.style.fontFamily = 'sans-serif';
  title.style.fontSize = '16px';
  title.style.fontWeight = '600';
  title.style.color = '#222';
  title.style.background = 'rgba(255,255,255,0.75)';
  title.style.padding = '4px 10px';
  title.style.borderRadius = '6px';
  title.style.pointerEvents = 'none';
  container.appendChild(title);

  var msgEl = panel.querySelector('[data-reportr=\"msg\"]');
  var inputEl = panel.querySelector('[data-reportr=\"regex\"]');
  var applyBtn = panel.querySelector('[data-reportr=\"apply\"]');
  var clearBtn = panel.querySelector('[data-reportr=\"clear\"]');
  var showLabelsEl = panel.querySelector('[data-reportr=\"showLabels\"]');
  var labelModeEl = panel.querySelector('[data-reportr=\"labelMode\"]');
  var samplesEl = panel.querySelector('[data-reportr=\"samples\"]');
  var tasksEl = panel.querySelector('[data-reportr=\"tasks\"]');
  var keepSamplesEl = panel.querySelector('[data-reportr=\"keepSamples\"]');
  var keepComparativesEl = panel.querySelector('[data-reportr=\"keepComparatives\"]');
  var minCovTypeEl = panel.querySelector('[data-reportr=\"minCovType\"]');
  var hideIsolatedEl = panel.querySelector('[data-reportr=\"hideIsolated\"]');
  var applyFiltersBtn = panel.querySelector('[data-reportr=\"applyFilters\"]');
  var resetFiltersBtn = panel.querySelector('[data-reportr=\"resetFilters\"]');
  var ttSizeEl = panel.querySelector('[data-reportr=\"tt_size\"]');
  var ttPdivEl = panel.querySelector('[data-reportr=\"tt_pdiv\"]');
  var ttTareanEl = panel.querySelector('[data-reportr=\"tt_tarean\"]');

  function setMsg(s) { if (msgEl) msgEl.textContent = s; }
  setMsg('Initializing graph…');

  function resolveNetwork() {
    // visNetwork binding stores the vis.Network instance on the inner div graph+el.id.
    // Depending on htmlwidgets versions, HTMLWidgets.find(...).chart may be absent.
    try {
      var graphEl = document.getElementById('graph' + el.id);
      if (graphEl && graphEl.chart) return graphEl.chart;
    } catch (e) {}

    var widget = null;
    try {
      widget = (window.HTMLWidgets && HTMLWidgets.find) ? HTMLWidgets.find('#' + el.id) : null;
    } catch (e) {}
    if (!widget) return null;

    // Some builds expose a helper.
    try {
      if (typeof widget.getNetwork === 'function') {
        var n0 = widget.getNetwork();
        if (n0) return n0;
      }
    } catch (e) {}

    // Fallback: chart may be the network or a function returning it.
    try {
      if (typeof widget.chart === 'function') {
        var n1 = widget.chart();
        if (n1) return n1;
      }
      if (widget.chart) return widget.chart;
    } catch (e) {}
    return null;
  }

  function resolveNodesDS(network) {
    return (network && network.body && network.body.data) ? network.body.data.nodes : null;
  }

  function startInit() {
    // Large graphs can take time to fully initialize. Be patient and keep retrying.
    Promise.all([
      _waitFor(resolveNetwork, 60000, 100)
    ]).then(function(res) {
      var network = res[0];
      var edgesDS = (network && network.body && network.body.data) ? network.body.data.edges : null;
      // Ensure the vis.js canvas uses the full available space.
      try {
        if (typeof network.setSize === 'function') network.setSize('100%', '100%');
        if (typeof network.redraw === 'function') network.redraw();
      } catch (e) {}
      return _waitFor(function(){ return resolveNodesDS(network); }, 60000, 100).then(function(nodesDS) {
      applyBtn.disabled = false;
      clearBtn.disabled = false;
      setMsg('');

      // Persist original node/edge styles once (for clean restore on Clear).
      var orig = window.__reportrGraphOrig || { nodes: {}, edges: {} };
      window.__reportrGraphOrig = orig;

      function snapshotOrig(nodesDS, edgesDS) {
        try {
          var nodes = nodesDS.get();
          for (var i = 0; i < nodes.length; i++) {
            var n = nodes[i];
            if (!n || n.id === undefined || n.id === null) continue;
            if (orig.nodes[n.id] === undefined) {
              orig.nodes[n.id] = {
                color: n.color,
                font: n.font,
                hidden: n.hidden
              };
            }
          }
        } catch (e) {}
        try {
          var edges = edgesDS.get();
          for (var j = 0; j < edges.length; j++) {
            var ed = edges[j];
            if (!ed || ed.id === undefined || ed.id === null) continue;
            if (orig.edges[ed.id] === undefined) {
              orig.edges[ed.id] = {
                color: ed.color,
                width: ed.width,
                hidden: ed.hidden
              };
            }
          }
        } catch (e) {}
      }

      function restoreOrig(nodesDS, edgesDS) {
        try {
          var nodeUpdates = [];
          for (var nid in orig.nodes) {
            if (!Object.prototype.hasOwnProperty.call(orig.nodes, nid)) continue;
            var o = orig.nodes[nid];
            nodeUpdates.push({ id: nid, color: o.color, font: o.font, hidden: o.hidden });
          }
          if (nodeUpdates.length) nodesDS.update(nodeUpdates);
        } catch (e) {}
        try {
          var edgeUpdates = [];
          for (var eid in orig.edges) {
            if (!Object.prototype.hasOwnProperty.call(orig.edges, eid)) continue;
            var eo = orig.edges[eid];
            edgeUpdates.push({ id: eid, color: eo.color, width: eo.width, hidden: eo.hidden });
          }
          if (edgeUpdates.length) edgesDS.update(edgeUpdates);
        } catch (e) {}
      }

      function applyDim(nodesDS, edgesDS, keepIds) {
        var keep = {};
        for (var i = 0; i < keepIds.length; i++) keep[keepIds[i]] = true;

        // Nodes: dim everything else, emphasize kept.
        var nodeUpdates = [];
        var nodes = nodesDS.get();
        for (var j = 0; j < nodes.length; j++) {
          var n = nodes[j];
          if (!n || n.id === undefined || n.id === null) continue;
          var inKeep = !!keep[n.id];
          if (inKeep) {
            nodeUpdates.push({
              id: n.id,
              hidden: false,
              color: { border: \"#1f77b4\", background: \"rgba(31,119,180,0.25)\" },
              font: { color: \"#111\" }
            });
          } else {
            nodeUpdates.push({
              id: n.id,
              hidden: false,
              color: { border: \"rgba(170,170,170,0.35)\", background: \"rgba(220,220,220,0.12)\" },
              font: { color: \"rgba(120,120,120,0.6)\" }
            });
          }
        }
        nodesDS.update(nodeUpdates);

        // Edges: emphasize if both endpoints kept, otherwise dim strongly.
        var edgeUpdates = [];
        var edges = edgesDS.get();
        for (var k = 0; k < edges.length; k++) {
          var ed = edges[k];
          if (!ed || ed.id === undefined || ed.id === null) continue;
          var fromKeep = !!keep[ed.from];
          var toKeep = !!keep[ed.to];
          var strong = fromKeep && toKeep;
          edgeUpdates.push({
            id: ed.id,
            hidden: false,
            width: strong ? 2.0 : 0.5,
            color: strong ? \"rgba(0,0,0,0.55)\" : \"rgba(180,180,180,0.18)\"
          });
        }
        edgesDS.update(edgeUpdates);
      }

      function _parseSampleId(nodeId) {
        var m = String(nodeId || '').match(/SMPL=([^|]+)/);
        return m ? m[1] : '';
      }

      function _covRank(ct) {
        ct = String(ct || '').toLowerCase();
        if (ct === 'near-full') return 4;
        if (ct === 'composite') return 3;
        if (ct === 'partial') return 2;
        if (ct === 'weak') return 1;
        return 0;
      }

      function _buildTitle(n) {
        // Base fields always present.
        var out = '';
        out += '<b>' + String(n.label_repeat || n.id || '') + '</b><br/>';
        if (n.cluster_id) out += 'cluster: ' + String(n.cluster_id) + '<br/>';
        if (n.subject_type || n.subject_id) out += 'subject: ' + String(n.subject_type || '') + ' ' + String(n.subject_id || '') + '<br/>';

        if (ttSizeEl && ttSizeEl.checked && n.size_pct !== undefined && n.size_pct !== null && String(n.size_pct) !== '') {
          out += 'size_pct: ' + String(n.size_pct) + '<br/>';
        }
        if (ttPdivEl && ttPdivEl.checked && n.pdiv_mean !== undefined && n.pdiv_mean !== null && String(n.pdiv_mean) !== '') {
          out += 'pdiv_mean: ' + String(n.pdiv_mean) + '<br/>';
        }
        if (ttTareanEl && ttTareanEl.checked && n.TAREAN_annotation) {
          out += 'TAREAN: ' + String(n.TAREAN_annotation) + '<br/>';
        }

        if (n.best_hit) out += 'NCBI: ' + String(n.best_hit) + (n.best_species ? (' (' + String(n.best_species) + ')') : '') + '<br/>';
        if (n.best_title) out += 'NCBI desc: ' + String(n.best_title) + '<br/>';
        if (n.oligo_best_sseqid) out += 'Oligo: ' + String(n.oligo_best_sseqid) + (n.oligo_fitting ? (' (' + String(n.oligo_fitting) + ')') : '') + '<br/>';
        return out;
      }

      function _applyLabelMode(nodesDS, showLabels, mode) {
        var nodes = nodesDS.get();
        var upd = [];
        for (var i = 0; i < nodes.length; i++) {
          var n = nodes[i];
          if (!n || n.id === undefined || n.id === null) continue;
          var label = '';
          if (showLabels) {
            var rep = String(n.label_repeat || n.id || '');
            var size = (n.size_pct !== undefined && n.size_pct !== null) ? String(n.size_pct) : '';
            var pdiv = (n.pdiv_mean !== undefined && n.pdiv_mean !== null) ? String(n.pdiv_mean) : '';
            if (mode === 'repeat') label = rep;
            else if (mode === 'repeat_size') label = size ? (rep + ' [' + size + '%]') : rep;
            else if (mode === 'repeat_pdiv') label = pdiv ? (rep + ' [pdiv=' + pdiv + ']') : rep;
            else if (mode === 'size_only') label = size ? (size + '%') : '';
            else label = rep;
          }
          upd.push({ id: n.id, label: label });
        }
        nodesDS.update(upd);
      }

      function _refreshTitles(nodesDS) {
        var nodes = nodesDS.get();
        var upd = [];
        for (var i = 0; i < nodes.length; i++) {
          var n = nodes[i];
          if (!n || n.id === undefined || n.id === null) continue;
          upd.push({ id: n.id, title: _buildTitle(n) });
        }
        nodesDS.update(upd);
      }

      function _populateSelect(selectEl, values) {
        if (!selectEl) return;
        selectEl.innerHTML = '';
        for (var i = 0; i < values.length; i++) {
          var opt = document.createElement('option');
          opt.value = values[i];
          opt.textContent = values[i];
          opt.selected = true;
          selectEl.appendChild(opt);
        }
      }

      function _getSelected(selectEl) {
        var out = [];
        if (!selectEl) return out;
        for (var i = 0; i < selectEl.options.length; i++) {
          var o = selectEl.options[i];
          if (o.selected) out.push(o.value);
        }
        return out;
      }

      function applyRegex() {
        var raw = (inputEl && inputEl.value) ? inputEl.value : '';
        raw = String(raw || '').trim();
        if (!raw) return setMsg('Enter a regex.');
        var re;
        try { re = _compileRegex(raw); } catch (e) { return setMsg('Invalid regex: ' + String(e)); }

        var nodes = nodesDS.get();
        var matched = _findMatchingNodeIds(nodes, re);
        if (matched.length === 0) return setMsg('No matches.');

        var ids = _expandWithNeighbors(network, matched);
        network.selectNodes(ids, false);
        if (edgesDS) {
          snapshotOrig(nodesDS, edgesDS);
          applyDim(nodesDS, edgesDS, ids);
        }
        if (matched.length === 1) network.focus(matched[0], { scale: 1.2, animation: { duration: 300 } });
        else network.fit({ nodes: ids, animation: { duration: 300 } });
        setMsg('Matched ' + matched.length + ' node(s); selected ' + ids.length + ' incl. neighbors.');
      }

      applyBtn.addEventListener('click', applyRegex);
      inputEl.addEventListener('keydown', function(e) { if (e.key === 'Enter') applyRegex(); });
      clearBtn.addEventListener('click', function() {
        inputEl.value = '';
        setMsg('');
        if (edgesDS) restoreOrig(nodesDS, edgesDS);
        network.unselectAll();
        network.fit({ animation: { duration: 200 } });
      });

      // Initialize control lists from current data.
      try {
        var allNodes = nodesDS.get();
        var smpls = {};
        for (var i0 = 0; i0 < allNodes.length; i0++) {
          var s0 = _parseSampleId(allNodes[i0].id);
          if (s0) smpls[s0] = true;
        }
        _populateSelect(samplesEl, Object.keys(smpls).sort());
      } catch (e) {}

      try {
        if (edgesDS) {
          var allEdges = edgesDS.get();
          var tasks = {};
          for (var j0 = 0; j0 < allEdges.length; j0++) {
            var t0 = String(allEdges[j0].task || '').trim();
            if (t0) tasks[t0] = true;
          }
          _populateSelect(tasksEl, Object.keys(tasks).sort());
        }
      } catch (e) {}

      function _applyFilters() {
        var selectedSamples = _getSelected(samplesEl);
        var selectedTasks = _getSelected(tasksEl);
        var keepSamples = keepSamplesEl ? !!keepSamplesEl.checked : true;
        var keepComparatives = keepComparativesEl ? !!keepComparativesEl.checked : true;
        var minCov = minCovTypeEl ? String(minCovTypeEl.value || 'weak') : 'weak';
        var minCovRank = _covRank(minCov);

        var keepSampleSet = {};
        for (var i = 0; i < selectedSamples.length; i++) keepSampleSet[selectedSamples[i]] = true;
        var keepTaskSet = {};
        for (var k = 0; k < selectedTasks.length; k++) keepTaskSet[selectedTasks[k]] = true;

        var nodeArr = nodesDS.get();
        var keepNode = {};
        var nodeUpdates = [];
        for (var n = 0; n < nodeArr.length; n++) {
          var nd = nodeArr[n];
          var isComp = String(nd.subject_type || '') === 'comparative';
          var keep = true;
          if (isComp && !keepComparatives) keep = false;
          if (!isComp && keepSamples) {
            var sid = _parseSampleId(nd.id);
            if (sid && !keepSampleSet[sid]) keep = false;
          }
          keepNode[nd.id] = !!keep;
          nodeUpdates.push({ id: nd.id, hidden: !keep });
        }
        nodesDS.update(nodeUpdates);

        // Edge filtering and (optional) isolation trimming.
        var edgeArr = edgesDS ? edgesDS.get() : [];
        var edgeUpdates = [];
        var degree = {};
        for (var e = 0; e < edgeArr.length; e++) {
          var ed = edgeArr[e];
          var from = ed.from, to = ed.to;
          var ok = true;
          if (!keepNode[from] || !keepNode[to]) ok = false;
          var task = String(ed.task || '').trim();
          if (ok && selectedTasks.length && !keepTaskSet[task]) ok = false;
          var ct = String(ed.coverage_type || '').trim();
          if (ok && _covRank(ct) < minCovRank) ok = false;
          edgeUpdates.push({ id: ed.id, hidden: !ok });
          if (ok) {
            degree[from] = (degree[from] || 0) + 1;
            degree[to] = (degree[to] || 0) + 1;
          }
        }
        if (edgesDS) edgesDS.update(edgeUpdates);

        if (hideIsolatedEl && hideIsolatedEl.checked) {
          var nu = [];
          for (var n2 = 0; n2 < nodeArr.length; n2++) {
            var nd2 = nodeArr[n2];
            var keep2 = keepNode[nd2.id] && (degree[nd2.id] || 0) > 0;
            nu.push({ id: nd2.id, hidden: !keep2 });
            keepNode[nd2.id] = keep2;
          }
          nodesDS.update(nu);
        }

        _applyLabelMode(nodesDS, showLabelsEl && showLabelsEl.checked, labelModeEl ? labelModeEl.value : 'repeat');
        _refreshTitles(nodesDS);
        setMsg('Filters applied.');
      }

      function _resetFilters() {
        if (samplesEl) {
          for (var i = 0; i < samplesEl.options.length; i++) samplesEl.options[i].selected = true;
        }
        if (tasksEl) {
          for (var j = 0; j < tasksEl.options.length; j++) tasksEl.options[j].selected = true;
        }
        if (keepSamplesEl) keepSamplesEl.checked = true;
        if (keepComparativesEl) keepComparativesEl.checked = true;
        if (minCovTypeEl) minCovTypeEl.value = 'near-full';
        if (hideIsolatedEl) hideIsolatedEl.checked = false;
        if (showLabelsEl) showLabelsEl.checked = false;
        if (labelModeEl) labelModeEl.value = 'repeat';

        // Unhide everything, then reapply title/label defaults.
        var nodesAll = nodesDS.get();
        var nu = [];
        for (var n = 0; n < nodesAll.length; n++) nu.push({ id: nodesAll[n].id, hidden: false, label: '' });
        nodesDS.update(nu);
        if (edgesDS) {
          var edgesAll = edgesDS.get();
          var eu = [];
          for (var e = 0; e < edgesAll.length; e++) eu.push({ id: edgesAll[e].id, hidden: false });
          edgesDS.update(eu);
        }
        _refreshTitles(nodesDS);
        setMsg('Filters reset.');
      }

      // Enable filter buttons now that datasets exist.
      applyFiltersBtn.disabled = false;
      resetFiltersBtn.disabled = false;

      applyFiltersBtn.addEventListener('click', _applyFilters);
      resetFiltersBtn.addEventListener('click', _resetFilters);

      // Live refresh of titles when tooltip field toggles change.
      if (ttSizeEl) ttSizeEl.addEventListener('change', function(){ _refreshTitles(nodesDS); });
      if (ttPdivEl) ttPdivEl.addEventListener('change', function(){ _refreshTitles(nodesDS); });
      if (ttTareanEl) ttTareanEl.addEventListener('change', function(){ _refreshTitles(nodesDS); });
      if (showLabelsEl) showLabelsEl.addEventListener('change', function(){ _applyLabelMode(nodesDS, showLabelsEl.checked, labelModeEl ? labelModeEl.value : 'repeat'); });
      if (labelModeEl) labelModeEl.addEventListener('change', function(){ _applyLabelMode(nodesDS, showLabelsEl && showLabelsEl.checked, labelModeEl.value); });

      // Initial tooltip content.
      _refreshTitles(nodesDS);
      });
    }).catch(function() {
      setMsg('Graph not ready yet (network not initialized). Waiting…');
      setTimeout(startInit, 500);
    });
  }
  startInit();
}
")

saveWidget(g, file = out_html, selfcontained = FALSE, title = "RepOrtR graph")

