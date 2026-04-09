/**
 * Celestial Codon — Frontend JavaScript
 * =======================================
 * Handles all UI interactions, API calls, and rendering for the DNA analysis platform.
 *
 * Key sections:
 *  - analyze()         : Sends DNA sequence to /analyze, renders all result tabs
 *  - renderAll()       : Orchestrates rendering of every tab after analysis
 *  - renderOverview()  : Stats, sequence viewer, restriction sites, palindromes
 *  - renderORFs()      : 6-frame ORF canvas map + ORF list
 *  - renderProtein()   : Protein stats, domain map, amino acid composition
 *  - renderCodons()    : Codon bias comparison table
 *  - renderMotifLib()  : Motif library scan results
 *  - renderCharts()    : Matplotlib chart images (nucleotide freq, GC window, hydrophobicity)
 *  - renderAIInsights(): Rule-based AI insight cards (no API needed)
 *  - smartBioAI()      : Local logic AI — answers common bio questions instantly
 *  - askGroq()         : Calls /ask_groq backend for complex questions (Groq free tier)
 *  - sendAiMessage()   : Chat panel handler with typewriter effect
 *  - exportXlsx/Fasta/Json(): Download result files
 *
 * AI Flow:
 *  User message → smartBioAI() [instant, no network]
 *              → if __NO_ANSWER__ → askGroq() [Groq API via backend]
 */

let currentSeq = '';
let currentData = null;

// ── Theme ──────────────────────────────────────────────────────────────────
function toggleTheme() {
  const html = document.documentElement;
  const isDark = html.getAttribute('data-theme') === 'dark';
  html.setAttribute('data-theme', isDark ? 'light' : 'dark');
  document.getElementById('themeBtn').textContent = isDark ? 'Light' : 'Dark';
}

// ── Advanced options toggle ────────────────────────────────────────────────
function toggleAdvanced() {
  const panel = document.getElementById('advPanel');
  const arrow = document.getElementById('advArrowSvg');
  const open = !panel.classList.contains('hidden');
  panel.classList.toggle('hidden', open);
  if (arrow) arrow.style.transform = open ? '' : 'rotate(90deg)';
}

// ── Sidebar nav ────────────────────────────────────────────────────────────
function showTab(name, btn) {
  document.querySelectorAll('.tab-content').forEach(el => el.classList.add('hidden'));
  document.querySelectorAll('.nav-item').forEach(el => el.classList.remove('active'));
  document.getElementById('tab-' + name).classList.remove('hidden');
  if (btn) btn.classList.add('active');
}

function showTabByName(name) {
  document.querySelectorAll('.tab-content').forEach(el => el.classList.add('hidden'));
  document.querySelectorAll('.nav-item').forEach(el => el.classList.remove('active'));
  document.getElementById('tab-' + name).classList.remove('hidden');
  document.querySelectorAll('.nav-item').forEach(btn => {
    if (btn.getAttribute('onclick') && btn.getAttribute('onclick').includes("'" + name + "'")) {
      btn.classList.add('active');
    }
  });
}

// ── Samples ────────────────────────────────────────────────────────────────
function loadSample(seq) { document.getElementById('seq').value = seq; }

// ── Analyze ────────────────────────────────────────────────────────────────
async function analyze() {
  console.log('analyze() triggered');

  const seq = document.getElementById('seq').value.trim();
  if (!seq) {
    showError('Please enter a DNA sequence before analyzing.');
    return;
  }

  const motif = document.getElementById('motif') ? document.getElementById('motif').value.trim() : '';
  const organism = document.getElementById('organism') ? document.getElementById('organism').value : 'ecoli';
  const theme = document.documentElement.getAttribute('data-theme');
  const btn = document.getElementById('analyzeBtn');

  document.getElementById('error').classList.add('hidden');
  btn.disabled = true;

  // Show loading overlay
  const overlay = document.getElementById('loadingOverlay');
  const stepEl = document.getElementById('loadingStep');
  const barEl = document.getElementById('loadingBar');

  if (overlay) {
    overlay.classList.add('visible');
    barEl.style.width = '10%';
  }

  const steps = [
    {text: 'Analyzing sequence...', pct: 25},
    {text: 'Detecting ORFs...', pct: 50},
    {text: 'Calculating protein properties...', pct: 75},
    {text: 'Generating insights...', pct: 90},
  ];
  let stepIdx = 0;
  const stepInterval = setInterval(() => {
    if (stepIdx < steps.length && stepEl && barEl) {
      stepEl.textContent = steps[stepIdx].text;
      barEl.style.width = steps[stepIdx].pct + '%';
      stepIdx++;
    }
  }, 600);

  const hideOverlay = () => {
    clearInterval(stepInterval);
    if (overlay) overlay.classList.remove('visible');
    if (barEl) barEl.style.width = '0%';
    btn.disabled = false;
  };

  try {
    console.log('Sending fetch to /analyze with sequence length:', seq.length);
    const res = await fetch('/analyze', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({sequence: seq, motif, organism, theme})
    });

    console.log('Response status:', res.status);
    const data = await res.json();
    console.log('Response data keys:', Object.keys(data));

    if (barEl) barEl.style.width = '100%';
    if (stepEl) stepEl.textContent = 'Complete';

    setTimeout(() => {
      hideOverlay();
      if (data.error) {
        showError(data.error);
        return;
      }
      currentSeq = data.sequence;
      currentData = data;
      console.log('Calling renderAll...');
      renderAll(data);
      document.getElementById('exportRow').style.display = '';
      loadSecondary();
    }, 300);

  } catch(e) {
    console.error('analyze() error:', e);
    hideOverlay();
    showError('Connection error: ' + e.message + '. Is the server running at http://127.0.0.1:5000?');
  }
}

function showError(msg) {
  const el = document.getElementById('error');
  el.textContent = '⚠ ' + msg;
  el.classList.remove('hidden');
  el.scrollIntoView({behavior: 'smooth', block: 'nearest'});
}

// ── Render all ─────────────────────────────────────────────────────────────
function renderAll(d) {
  console.log('renderAll() called');
  try { renderOverview(d); } catch(e) { console.error('renderOverview failed:', e); }
  try { renderORFs(d); } catch(e) { console.error('renderORFs failed:', e); }
  try { renderProtein(d); } catch(e) { console.error('renderProtein failed:', e); }
  try { renderCodons(d); } catch(e) { console.error('renderCodons failed:', e); }
  try { renderMotifLib(d); } catch(e) { console.error('renderMotifLib failed:', e); }
  try { renderCharts(d); } catch(e) { console.error('renderCharts failed:', e); }
  try { buildEnzymeCheckboxes(); } catch(e) { console.error('buildEnzymeCheckboxes failed:', e); }
  try { renderAIInsights(d); } catch(e) { console.error('renderAIInsights failed:', e); }

  const analysisNav = document.getElementById('analysisNav');
  const labNav = document.getElementById('labNav');
  const insightsNav = document.getElementById('insightsNav');
  if (analysisNav) analysisNav.style.display = '';
  if (labNav) labNav.style.display = '';
  if (insightsNav) insightsNav.style.display = '';

  // Show overview tab
  document.querySelectorAll('.tab-content').forEach(el => el.classList.add('hidden'));
  document.querySelectorAll('.nav-item').forEach(el => el.classList.remove('active'));
  const overviewTab = document.getElementById('tab-overview');
  if (overviewTab) overviewTab.classList.remove('hidden');
  document.querySelectorAll('.nav-item').forEach(btn => {
    if (btn.getAttribute('onclick') && btn.getAttribute('onclick').includes("'overview'")) {
      btn.classList.add('active');
    }
  });

  console.log('renderAll() complete — overview tab shown');
  // Update AI chat greeting if panel has been opened
  updateAiGreeting();
}

// ── AI Insights (rule-based, card layout) ─────────────────────────────────
function renderAIInsights(d) {
  const el = document.getElementById('r-aiInsights');
  if (!el) return;

  const gc = parseFloat(d.gc);
  const tm = parseFloat(d.tm);
  const totalOrfs = d.orfs_6frame ? d.orfs_6frame.reduce((s,f) => s + f.orfs.length, 0) : 0;
  const cards = [];
  const nextSteps = [];

  // GC content card
  if (gc > 60) {
    cards.push({type:'success', icon:'G', text:`GC content is high (${gc}%)`, detail:'Indicates strong thermal stability and resistance to denaturation. Suitable for high-temperature PCR applications.'});
  } else if (gc < 40) {
    cards.push({type:'warning', icon:'G', text:`GC content is low (${gc}%)`, detail:'Sequence may be AT-rich, common in regulatory regions or repeat elements. Lower melting temperature expected.'});
  } else {
    cards.push({type:'info', icon:'G', text:`GC content is moderate (${gc}%)`, detail:'Balanced sequence composition with typical thermal stability for standard molecular biology applications.'});
  }

  // Melting temp card
  if (tm > 80) {
    cards.push({type:'warning', icon:'T', text:`High melting temperature (${tm}°C)`, detail:'Primers designed from this sequence will require elevated annealing temperatures (typically Tm - 5°C).'});
  } else if (tm < 50) {
    cards.push({type:'warning', icon:'T', text:`Low melting temperature (${tm}°C)`, detail:'Consider using longer primers or adding DMSO for PCR optimization to improve specificity.'});
  } else {
    cards.push({type:'info', icon:'T', text:`Melting temperature is ${tm}°C`, detail:'Standard annealing temperature range. Suitable for conventional PCR protocols.'});
  }

  // ORF card
  if (totalOrfs > 0) {
    cards.push({type:'success', icon:'O', text:`${totalOrfs} open reading frame(s) detected`, detail:'Sequence contains protein-coding potential. Review the 6-Frame ORFs tab to identify the longest coding region.'});
    nextSteps.push('Review 6-Frame ORFs to identify the longest coding region');
  } else {
    cards.push({type:'danger', icon:'O', text:'No ORFs detected', detail:'Sequence may be non-coding, regulatory, or too short for translation. Try a longer sequence.'});
  }

  // Protein card
  if (d.protein_props) {
    const p = d.protein_props;
    const hydro = parseFloat(p.avg_hydrophobicity);
    if (hydro > 0.5) {
      cards.push({type:'warning', icon:'P', text:`Protein is hydrophobic (avg ${hydro})`, detail:'May interact with lipid membranes or form transmembrane domains. Consider detergent-based extraction.'});
    } else if (hydro < -1) {
      cards.push({type:'success', icon:'P', text:`Protein is hydrophilic (avg ${hydro})`, detail:'Likely soluble and cytoplasmic. Good candidate for recombinant expression in E. coli.'});
    } else {
      cards.push({type:'info', icon:'P', text:`Protein has mixed hydrophobicity (avg ${hydro})`, detail:`MW: ${p.mw.toLocaleString()} Da, pI: ${p.pi}. Check the Protein tab for full domain analysis.`});
    }
    nextSteps.push(`Protein MW is ${p.mw.toLocaleString()} Da with pI ${p.pi} — check Protein tab`);
  }

  // Restriction sites
  const rsites = Object.keys(d.restriction_sites || {});
  if (rsites.length > 0) {
    cards.push({type:'info', icon:'R', text:`Restriction sites: ${rsites.join(', ')}`, detail:'These sites are useful for cloning strategy design. Use Lab Tools to simulate enzyme digestion.'});
    nextSteps.push('Use Lab Tools to simulate restriction digest');
  }

  // Motifs
  if (d.motif_library && d.motif_library.length > 0) {
    cards.push({type:'success', icon:'M', text:`${d.motif_library.length} regulatory motif(s) detected`, detail:'Sequence may contain promoter or transcription factor binding elements. Check the Motifs tab for details.'});
  }

  nextSteps.push('Run the Mutator to simulate point mutations and check disease associations');
  nextSteps.push('Export results as Excel for detailed offline analysis');

  let cardIdx = 0;
  el.innerHTML = `<div class="ai-cards">` +
    cards.map(c => {
      const idx = cardIdx++;
      return `<div class="ai-card ${c.type}">
        <div class="ai-card-icon">${c.icon}</div>
        <div class="ai-card-body">
          <div class="ai-card-text">${c.text}</div>
          <div class="ai-card-detail" id="ai-detail-${idx}">${c.detail}</div>
          <button class="ai-card-toggle" onclick="toggleAiDetail(${idx}, this)">Show more</button>
        </div>
      </div>`;
    }).join('') +
    `</div>
    <div class="ai-next-steps">
      <div class="ai-next-title">Suggested Next Steps</div>
      ${nextSteps.map(s => `<div class="ai-next-item">${s}</div>`).join('')}
    </div>`;

  // Also populate the preview on overview
  const previewEl = document.getElementById('r-aiPreview');
  if (previewEl) {
    previewEl.innerHTML = `<div class="ai-cards">` +
      cards.slice(0, 3).map(c => `<div class="ai-card ${c.type}" style="padding:.65rem .85rem">
        <div class="ai-card-icon" style="width:22px;height:22px;font-size:.72rem">${c.icon}</div>
        <div class="ai-card-text" style="font-size:.8rem">${c.text}</div>
      </div>`).join('') +
      `</div>`;
  }
}

function toggleAiDetail(idx, btn) {
  const detail = document.getElementById('ai-detail-' + idx);
  const open = detail.classList.toggle('open');
  btn.textContent = open ? 'Show less' : 'Show more';
}

async function handleAskAi() {
  const input = document.getElementById('askAiInput');
  const resp = document.getElementById('r-askAiResponse');
  const q = input.value.trim();
  if (!q) return;
  input.disabled = true;
  resp.textContent = 'Thinking...';
  resp.style.color = 'var(--text-muted)';
  try {
    const answer = await askAI(q);
    resp.textContent = answer || 'No response.';
    resp.style.color = 'var(--text)';
  } catch (err) {
    resp.textContent = 'Error: ' + err.message;
    resp.style.color = 'var(--red)';
  } finally {
    input.disabled = false;
  }
}

// ── Overview ───────────────────────────────────────────────────────────────
function renderOverview(d) {
  document.getElementById('r-length').textContent = d.length;
  document.getElementById('r-gc').textContent = d.gc + '%';
  document.getElementById('r-at').textContent = d.at + '%';
  document.getElementById('r-tm').textContent = d.tm + '°C';
  document.getElementById('overview-subtitle').textContent =
    `Sequence: ${d.sequence.slice(0, 40)}${d.sequence.length > 40 ? '…' : ''}`;

  // Sequence viewer
  const colors = {A:'base-A', T:'base-T', C:'base-C', G:'base-G'};
  let html = '';
  for (let i = 0; i < d.sequence.length; i++) {
    if (i % 10 === 0) html += `<span class="seq-pos">${i+1}</span>`;
    const b = d.sequence[i];
    html += `<span class="${colors[b]||''}" title="pos ${i+1}">${b}</span>`;
    if ((i+1) % 60 === 0) html += '<br>';
  }
  document.getElementById('r-seqviewer').innerHTML = html;

  // Freq bars
  const bColors = {A:'#3fb950', T:'#f85149', C:'#58a6ff', G:'#e3b341'};
  const max = Math.max(...Object.values(d.freq));
  document.getElementById('r-freq').innerHTML = Object.entries(d.freq).map(([b, n]) => `
    <div class="freq-item">
      <span class="freq-base" style="color:${bColors[b]}">${b}</span>
      <div class="freq-track"><div class="freq-fill" style="width:${n/max*100}%;background:${bColors[b]}"></div></div>
      <span class="freq-count">${n}</span>
    </div>`).join('');

  document.getElementById('r-rc').textContent = d.reverse_complement;
  document.getElementById('r-rna').textContent = d.rna;
  document.getElementById('r-protein').textContent = d.protein;

  // Motif
  const mc = document.getElementById('motif-card');
  if (d.motif) {
    mc.style.display = '';
    const pos = d.motif_positions;
    document.getElementById('r-motif').innerHTML = pos.length
      ? `<p>Found <strong>${pos.length}</strong> occurrence(s) of <code style="display:inline;padding:.1rem .4rem">${d.motif}</code> at positions: <span class="positions">${pos.join(', ')}</span></p>`
      : `<p class="none">Motif "${d.motif}" not found in sequence.</p>`;
  } else { mc.style.display = 'none'; }

  // Restriction sites — table
  const sites = Object.entries(d.restriction_sites);
  document.getElementById('r-restriction').innerHTML = sites.length
    ? `<table class="rest-table"><thead><tr><th>Enzyme</th><th>Site</th><th>Positions</th></tr></thead><tbody>
        ${sites.map(([e, i]) => `<tr>
          <td><span class="enzyme-name">${e}</span></td>
          <td><span class="site-badge">${i.site}</span></td>
          <td class="positions">${i.positions.join(', ')}</td>
        </tr>`).join('')}
      </tbody></table>`
    : '<p class="none">No common restriction sites found.</p>';

  // Palindromes
  document.getElementById('r-palindromes').innerHTML = d.palindromes.length
    ? d.palindromes.map(p => `<span class="pal-tag" title="pos ${p.start}–${p.end}">${p.seq}</span>`).join('')
    : '<p class="none">No palindromic sequences found.</p>';
}

// ── 6-Frame ORFs ───────────────────────────────────────────────────────────
function renderORFs(d) {
  const seqLen = d.length;
  const canvas = document.getElementById('orfCanvas');
  const W = canvas.parentElement.clientWidth || 760;
  canvas.width = W; canvas.height = 210;
  const ctx = canvas.getContext('2d');
  const isDark = document.documentElement.getAttribute('data-theme') === 'dark';
  ctx.fillStyle = isDark ? '#060d1a' : '#f6f8fa';
  ctx.fillRect(0, 0, W, 210);

  const frameColors = ['#388bfd','#bc8cff','#3fb950','#f85149','#f0883e','#e3b341'];
  const frameLabels = ['+1','+2','+3','-1','-2','-3'];
  const rowH = 22, rowPad = 8, topPad = 16, labelW = 28;

  d.orfs_6frame.forEach((frame, fi) => {
    const y = topPad + fi * (rowH + rowPad);
    ctx.fillStyle = isDark ? '#0f1e35' : '#e8eef5';
    ctx.fillRect(labelW + 4, y, W - labelW - 8, rowH);
    ctx.fillStyle = frameColors[fi];
    ctx.font = 'bold 10px monospace';
    ctx.fillText(frameLabels[fi], 2, y + rowH / 2 + 4);
    frame.orfs.forEach(orf => {
      const x = labelW + 4 + (orf.start / seqLen) * (W - labelW - 8);
      const w = Math.max(4, (orf.length / seqLen) * (W - labelW - 8));
      ctx.fillStyle = frameColors[fi];
      ctx.fillRect(x, y + 2, w, rowH - 4);
    });
  });

  let html = '';
  d.orfs_6frame.forEach((frame, fi) => {
    if (!frame.orfs.length) return;
    html += `<div class="frame-group">
      <span class="frame-label" style="color:${frameColors[fi]}">${frame.strand} Frame ${frame.frame}</span>`;
    frame.orfs.forEach(orf => {
      html += `<div class="orf-row">
        <span>Start: <strong>${orf.start}</strong></span>
        <span>End: <strong>${orf.end}</strong></span>
        <span>Length: <strong>${orf.length} bp</strong></span>
        ${orf.is_longest ? '<span class="badge-longest">Longest ORF</span>' : ''}
        <span class="orf-seq">${orf.seq.length > 80 ? orf.seq.slice(0,80)+'…' : orf.seq}</span>
      </div>`;
    });
    html += '</div>';
  });
  document.getElementById('r-orfs6').innerHTML = html || '<p class="none">No ORFs found.</p>';
}

// ── Protein ────────────────────────────────────────────────────────────────
function renderProtein(d) {
  const p = d.protein_props;
  const statsEl = document.getElementById('r-protStats');
  if (!p) { statsEl.innerHTML = '<p class="none">No protein translated.</p>'; return; }

  statsEl.innerHTML = `
    <div class="stat-card s-blue"><div class="stat-val c-blue">${p.length}</div><div class="stat-label">Amino Acids</div></div>
    <div class="stat-card s-green"><div class="stat-val c-green">${p.mw.toLocaleString()}</div><div class="stat-label">MW (Da)</div></div>
    <div class="stat-card s-orange"><div class="stat-val c-orange">${p.pi}</div><div class="stat-label">Isoelectric Point (pI)</div></div>
    <div class="stat-card s-purple"><div class="stat-val c-purple">${p.avg_hydrophobicity}</div><div class="stat-label">Avg Hydrophobicity</div></div>`;

  const AA_HYDRO = {Ile:4.5,Val:4.2,Leu:3.8,Phe:2.8,Cys:2.5,Met:1.9,Ala:1.8,Gly:-0.4,Thr:-0.7,Ser:-0.8,Trp:-0.9,Tyr:-1.3,Pro:-1.6,His:-3.2,Glu:-3.5,Gln:-3.5,Asp:-3.5,Asn:-3.5,Lys:-3.9,Arg:-4.5};
  const protein = d.protein_list || [];
  document.getElementById('r-domainMap').innerHTML = protein.length
    ? protein.map(aa => {
        const h = AA_HYDRO[aa] || 0;
        return `<span class="dom-aa ${h > 0 ? 'dom-hydro' : 'dom-hydrophil'}" title="${aa} (${h})">${aa.slice(0,1)}</span>`;
      }).join('')
    : '<p class="none">No protein.</p>';

  document.getElementById('r-aa').innerHTML = (d.aa_composition || []).map(({aa, count}) =>
    `<div class="aa-item"><span class="aa-name">${aa}</span><span class="aa-count">${count}</span></div>`
  ).join('') || '<p class="none">No amino acids.</p>';
}

// ── Codon Bias ─────────────────────────────────────────────────────────────
function renderCodons(d) {
  const tbody = document.getElementById('r-codons');
  if (!d.codon_bias || !d.codon_bias.length) {
    tbody.innerHTML = '<tr><td colspan="6" class="none" style="padding:.75rem">No codon data.</td></tr>'; return;
  }
  const maxDiff = Math.max(...d.codon_bias.map(c => Math.abs(c.diff)), 1);
  tbody.innerHTML = d.codon_bias.map(({codon, aa, seq_freq, ref_freq, diff}) => {
    const color = diff > 0 ? '#3fb950' : '#f85149';
    const barW = Math.abs(diff) / maxDiff * 100;
    return `<tr>
      <td><span class="codon-mono">${codon}</span></td>
      <td><span class="codon-aa-cell">${aa}</span></td>
      <td><span class="codon-num">${seq_freq}</span></td>
      <td><span class="codon-num">${ref_freq}</span></td>
      <td class="diff-bar-cell"><div class="diff-bar-wrap"><div class="diff-bar-fill" style="width:${barW}%;background:${color}"></div></div></td>
      <td><span class="diff-num" style="color:${color}">${diff > 0 ? '+' : ''}${diff}</span></td>
    </tr>`;
  }).join('');
}

// ── Motif Library ──────────────────────────────────────────────────────────
function renderMotifLib(d) {
  const el = document.getElementById('r-motiflib');
  if (!d.motif_library || !d.motif_library.length) {
    el.innerHTML = `<div class="empty-state">
      <svg class="empty-state-svg" width="56" height="56" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
        <circle cx="11" cy="11" r="8"/><line x1="21" y1="21" x2="16.65" y2="16.65"/>
        <line x1="11" y1="8" x2="11" y2="14"/><line x1="8" y1="11" x2="14" y2="11"/>
      </svg>
      <div class="empty-state-title">No motifs detected</div>
      <ul class="empty-state-tips">
        <li>Try using a longer sequence (200+ bp)</li>
        <li>Sequences with TATA-box, Kozak, or Shine-Dalgarno regions will show results</li>
        <li>The GFP sample sequence contains detectable motifs</li>
      </ul>
      <button class="empty-state-btn" onclick="loadSample('ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA');showTabByName('input')">Load GFP Example</button>
    </div>`; return;
  }
  el.innerHTML = d.motif_library.map(m => `
    <div class="motif-row">
      <span class="motif-name">${m.name}</span>
      <span class="site-badge">${m.pattern}</span>
      <span class="positions">${m.count} hit(s) at: ${m.positions.join(', ')}</span>
    </div>`).join('');
}

// ── Charts ─────────────────────────────────────────────────────────────────
function renderCharts(d) {
  document.getElementById('r-chart').src = 'data:image/png;base64,' + d.chart;
  const cc = document.getElementById('circularCard');
  if (d.circular) {
    cc.style.display = '';
    document.getElementById('r-circular').src = 'data:image/png;base64,' + d.circular;
  } else { cc.style.display = 'none'; }
}

// ── Lab Tools ──────────────────────────────────────────────────────────────
function buildEnzymeCheckboxes() {
  const enzymes = ['EcoRI','BamHI','HindIII','NotI','XhoI','NcoI','SalI','PstI','KpnI','SmaI','XbaI','SphI'];
  document.getElementById('enzymeCheckboxes').innerHTML = enzymes.map(e =>
    `<label class="enz-check"><input type="checkbox" value="${e}"> ${e}</label>`
  ).join('');

  if (currentData && currentData.primers) {
    const p = currentData.primers;
    document.getElementById('r-primers').innerHTML = `
      <div class="primer-grid">
        <div class="primer-card">
          <div class="primer-label">Forward Primer (5'→3')</div>
          <code>${p.forward.seq}</code>
          <div class="primer-stats">
            <span>Tm: <strong>${p.forward.tm}°C</strong></span>
            <span>GC: <strong>${p.forward.gc}%</strong></span>
            <span class="${p.forward.hairpin ? 'warn' : 'ok'}">${p.forward.hairpin ? '⚠ Hairpin risk' : '✓ No hairpin'}</span>
          </div>
        </div>
        <div class="primer-card">
          <div class="primer-label">Reverse Primer (5'→3')</div>
          <code>${p.reverse.seq}</code>
          <div class="primer-stats">
            <span>Tm: <strong>${p.reverse.tm}°C</strong></span>
            <span>GC: <strong>${p.reverse.gc}%</strong></span>
            <span class="${p.reverse.hairpin ? 'warn' : 'ok'}">${p.reverse.hairpin ? '⚠ Hairpin risk' : '✓ No hairpin'}</span>
          </div>
        </div>
      </div>
      <p style="margin-top:.75rem;color:var(--text-muted);font-size:.83rem">Expected amplicon: <strong>${p.amplicon_size} bp</strong></p>`;
  }
}

async function runDigest() {
  const checked = [...document.querySelectorAll('#enzymeCheckboxes input:checked')].map(e => e.value);
  if (!checked.length) { alert('Select at least one enzyme.'); return; }
  const res = await fetch('/digest', {method:'POST', headers:{'Content-Type':'application/json'},
    body: JSON.stringify({sequence: currentSeq, enzymes: checked})});
  const data = await res.json();
  const el = document.getElementById('r-digest');
  if (data.error || data.note) { el.innerHTML = `<p class="none">${data.error || data.note}</p>`; return; }
  el.innerHTML = `<p style="color:var(--text-muted);margin-bottom:.75rem;font-size:.83rem">Cut positions: <span class="positions">${data.cuts.join(', ')}</span></p>
    <div class="frag-list">${data.fragments.map((f, i) => `
      <div class="frag-row">
        <span class="frag-num">${i+1}</span>
        <span>${f.start}–${f.end}</span>
        <strong>${f.size} bp</strong>
      </div>`).join('')}</div>`;
}

async function runPCR() {
  const fwd = document.getElementById('pcrFwd').value.trim();
  const rev = document.getElementById('pcrRev').value.trim();
  const res = await fetch('/pcr', {method:'POST', headers:{'Content-Type':'application/json'},
    body: JSON.stringify({sequence: currentSeq, fwd, rev})});
  const data = await res.json();
  const el = document.getElementById('r-pcr');
  if (data.error) { el.innerHTML = `<p class="none">${data.error}</p>`; return; }
  if (!data.found) { el.innerHTML = `<p class="none">${data.message}</p>`; return; }
  el.innerHTML = `
    <div class="stats-row" style="grid-template-columns:repeat(3,1fr);margin-bottom:.75rem">
      <div class="stat-card s-blue"><div class="stat-val c-blue">${data.amplicon_size}</div><div class="stat-label">Amplicon (bp)</div></div>
      <div class="stat-card s-green"><div class="stat-val c-green">${data.amplicon_gc}%</div><div class="stat-label">GC Content</div></div>
      <div class="stat-card s-orange"><div class="stat-val c-orange">${data.fwd_pos}–${data.rev_pos}</div><div class="stat-label">Position</div></div>
    </div>
    <div class="seq-block"><div class="seq-block-label">Amplicon Sequence</div><code>${data.amplicon_seq}</code></div>`;
}

// ── Mutation ───────────────────────────────────────────────────────────────
async function runMutation() {
  const pos = parseInt(document.getElementById('mutPos').value);
  const base = document.getElementById('mutBase').value;
  if (!pos) { alert('Enter a position.'); return; }
  const res = await fetch('/mutate', {method:'POST', headers:{'Content-Type':'application/json'},
    body: JSON.stringify({sequence: currentSeq, position: pos, base})});
  const data = await res.json();
  const el = document.getElementById('r-mutation');
  if (data.error) { el.innerHTML = `<p class="none">${data.error}</p>`; return; }

  el.innerHTML = `
    <div class="mut-result">
      <div class="impact-badge" style="background:${data.impact_color}15;border-color:${data.impact_color};color:${data.impact_color}">
        ${data.impact} Mutation
      </div>
      <div class="mut-details">
        <div class="mut-detail"><div class="label">Position</div><div class="value">${data.position}</div></div>
        <div class="mut-detail"><div class="label">Base Change</div><div class="value">${data.original_base} → ${data.new_base}</div></div>
        <div class="mut-detail"><div class="label">Codon</div><div class="value">${data.original_codon} → ${data.mutated_codon}</div></div>
        <div class="mut-detail"><div class="label">Amino Acid</div><div class="value">${data.original_aa} → ${data.mutated_aa}</div></div>
      </div>
      <div class="seq-block"><div class="seq-block-label">Original Protein</div><code>${data.original_protein || '—'}</code></div>
      <div class="seq-block"><div class="seq-block-label">Mutated Protein</div><code class="protein">${data.mutated_protein || '—'}</code></div>
    </div>`;

  // Sequence comparison with highlighted mutation position
  const compareEl = document.getElementById('r-mutSeqCompare');
  if (compareEl) {
    compareEl.classList.remove('hidden');
    const colors = {A:'base-A', T:'base-T', C:'base-C', G:'base-G'};
    const mutIdx = pos - 1;
    const mutatedSeq = currentSeq.slice(0, mutIdx) + base + currentSeq.slice(mutIdx + 1);

    function renderSeqHighlight(seq, highlightIdx) {
      let html = '';
      for (let i = 0; i < seq.length; i++) {
        if (i % 10 === 0) html += `<span class="seq-pos">${i+1}</span>`;
        const b = seq[i];
        if (i === highlightIdx) {
          html += `<span class="mut-highlight" title="pos ${i+1} — mutated">${b}</span>`;
        } else {
          html += `<span class="${colors[b]||''}" title="pos ${i+1}">${b}</span>`;
        }
        if ((i+1) % 60 === 0) html += '<br>';
      }
      return html;
    }
    document.getElementById('r-origSeq').innerHTML = renderSeqHighlight(currentSeq, mutIdx);
    document.getElementById('r-mutSeq').innerHTML = renderSeqHighlight(mutatedSeq, mutIdx);
  }

  document.getElementById('diseaseCard').style.display = '';
  document.getElementById('r-disease').innerHTML = '';
}

async function runDiseaseLookup() {
  const pos = document.getElementById('mutPos').value;
  const base = document.getElementById('mutBase').value;
  const btn = document.getElementById('diseaseBtn');
  btn.disabled = true; btn.textContent = 'Querying database...';

  const res = await fetch('/disease', {method:'POST', headers:{'Content-Type':'application/json'},
    body: JSON.stringify({sequence: currentSeq, position: pos, base})});
  const data = await res.json();
  btn.disabled = false; btn.textContent = 'Query External Database';

  const el = document.getElementById('r-disease');
  if (data.error) { el.innerHTML = `<p class="none">${data.error}</p>`; return; }
  const d = data.disease_info;

  el.innerHTML = `
    <table class="disease-table">
      <tr><td class="dl">Mutation Type</td>
        <td><span class="impact-badge" style="background:${data.impact_color}15;border-color:${data.impact_color};color:${data.impact_color};padding:.2rem .65rem;font-size:.8rem">${data.impact}</span></td></tr>
      <tr><td class="dl">Codon Change</td>
        <td><code style="display:inline;padding:.1rem .4rem">${data.original_codon}</code> → <code style="display:inline;padding:.1rem .4rem">${data.mutated_codon}</code></td></tr>
      <tr><td class="dl">Amino Acid</td>
        <td><strong>${data.original_aa}</strong> → <strong>${data.mutated_aa}</strong></td></tr>
      ${d.variant_id ? `<tr><td class="dl">Variant ID</td><td style="font-family:'JetBrains Mono',monospace;font-size:.82rem;color:#79c0ff">${d.variant_id}</td></tr>` : ''}
      <tr><td class="dl">Disease</td>
        <td class="${d.found ? 'disease-found' : 'disease-none'}">${d.disease}</td></tr>
      <tr><td class="dl">Severity</td>
        <td><span class="sev-badge" style="background:${d.severity_color}15;border-color:${d.severity_color};color:${d.severity_color}">${d.severity}</span></td></tr>
      ${d.details.length ? `<tr><td class="dl">Predictions</td><td style="font-size:.82rem;color:var(--text-muted)">${d.details.join(' &nbsp;|&nbsp; ')}</td></tr>` : ''}
      ${d.clinvar_entries && d.clinvar_entries.length ? `<tr><td class="dl">ClinVar</td><td style="font-family:'JetBrains Mono',monospace;font-size:.8rem;color:#79c0ff">${d.clinvar_entries.join(' | ')}</td></tr>` : ''}
    </table>
    ${d.source ? `<p style="font-size:.72rem;color:var(--text-dim);margin-top:.5rem;text-align:right">Source: ${d.source}</p>` : ''}`;
}

// ── Secondary Structure ────────────────────────────────────────────────────
async function loadSecondary() {
  if (!currentSeq) return;
  const res = await fetch('/secondary', {method:'POST', headers:{'Content-Type':'application/json'},
    body: JSON.stringify({sequence: currentSeq})});
  const data = await res.json();
  if (data.error) return;

  const cpgEl = document.getElementById('r-cpg');
  cpgEl.innerHTML = data.cpg_islands.length
    ? data.cpg_islands.map(c => `<div class="orf-row">
        <span>Start: <strong>${c.start}</strong></span>
        <span>End: <strong>${c.end}</strong></span>
        <span>Length: <strong>${c.length} bp</strong></span>
      </div>`).join('')
    : '<p class="none">No CpG islands detected (sequence may be too short or GC-poor).</p>';

  const ssEl = document.getElementById('r-secondary');
  const ss = data.secondary_structure;
  if (!ss || !ss.length) { ssEl.innerHTML = '<p class="none">No protein to predict structure for.</p>'; return; }
  const labels = {H:'Alpha Helix', E:'Beta Sheet', C:'Coil'};
  ssEl.innerHTML = ss.map(r =>
    `<span class="ss-res ss-${r.struct}" title="${r.aa} (pos ${r.pos}) — ${labels[r.struct]}">${r.aa.slice(0,1)}</span>`
  ).join('');
}

// ── ELI5 ──────────────────────────────────────────────────────────────────
async function loadEli5() {
  const btn = document.getElementById('eli5Btn');
  btn.disabled = true; btn.textContent = 'Generating...';
  const res = await fetch('/eli5', {method:'POST', headers:{'Content-Type':'application/json'},
    body: JSON.stringify({sequence: currentSeq})});
  const data = await res.json();
  btn.disabled = false; btn.textContent = 'Regenerate';
  const el = document.getElementById('r-eli5');
  if (data.error) { el.innerHTML = `<p class="none">${data.error}</p>`; return; }
  el.innerHTML = data.lines.map((line, i) =>
    `<div class="eli5-item"><span class="eli5-num">${i+1}</span><p>${line}</p></div>`
  ).join('');
}

// ── Quiz ───────────────────────────────────────────────────────────────────
async function loadQuiz() {
  const btn = document.getElementById('quizBtn');
  btn.disabled = true; btn.textContent = 'Generating...';
  const res = await fetch('/quiz', {method:'POST', headers:{'Content-Type':'application/json'},
    body: JSON.stringify({sequence: currentSeq})});
  const data = await res.json();
  btn.disabled = false; btn.textContent = 'New Quiz';
  const el = document.getElementById('r-quiz');
  if (data.error) { el.innerHTML = `<p class="none">${data.error}</p>`; return; }
  el.innerHTML = data.questions.map((q, qi) => `
    <div class="quiz-item" id="quiz-${qi}">
      <p class="quiz-q"><strong>Q${qi+1}.</strong> ${q.q}</p>
      <div class="quiz-options">
        ${q.options.map((opt, oi) =>
          `<button class="quiz-opt" onclick="checkAnswer(${qi},${oi},${q.answer},'${q.explanation.replace(/'/g,"\\'")}')">
            <span class="quiz-letter">${'ABCD'[oi]}</span> ${opt}
          </button>`
        ).join('')}
      </div>
      <div class="quiz-feedback hidden" id="feedback-${qi}"></div>
    </div>`).join('');
}

function checkAnswer(qi, chosen, correct, explanation) {
  const card = document.getElementById(`quiz-${qi}`);
  const feedback = document.getElementById(`feedback-${qi}`);
  card.querySelectorAll('.quiz-opt').forEach((btn, i) => {
    btn.disabled = true;
    if (i === correct) btn.classList.add('correct');
    else if (i === chosen) btn.classList.add('wrong');
  });
  feedback.classList.remove('hidden');
  feedback.innerHTML = chosen === correct
    ? `<span class="ok">✓ Correct!</span> ${explanation}`
    : `<span class="warn">✗ Wrong.</span> ${explanation}`;
}

// ── Exports ────────────────────────────────────────────────────────────────
async function exportXlsx() {
  const res = await fetch('/export/xlsx', {method:'POST', headers:{'Content-Type':'application/json'},
    body: JSON.stringify({sequence: currentSeq})});
  const blob = await res.blob();
  const a = document.createElement('a'); a.href = URL.createObjectURL(blob);
  a.download = 'dna_report.xlsx'; a.click();
}

async function exportFasta() {
  const res = await fetch('/export/fasta', {method:'POST', headers:{'Content-Type':'application/json'},
    body: JSON.stringify({sequence: currentSeq, name: 'my_sequence'})});
  const blob = await res.blob();
  const a = document.createElement('a'); a.href = URL.createObjectURL(blob);
  a.download = 'sequence.fasta'; a.click();
}

async function exportJson() {
  const res = await fetch('/export/json', {method:'POST', headers:{'Content-Type':'application/json'},
    body: JSON.stringify({sequence: currentSeq})});
  const blob = await res.blob();
  const a = document.createElement('a'); a.href = URL.createObjectURL(blob);
  a.download = 'celestial_codon_report.json'; a.click();
}

// Ctrl+Enter shortcut
document.addEventListener('keydown', e => { if (e.ctrlKey && e.key === 'Enter') analyze(); });

// ── Floating AI Chat Panel ─────────────────────────────────────────────────
function toggleAiPanel() {
  const panel = document.getElementById('aiChatPanel');
  const isOpen = panel.classList.toggle('open');
  if (isOpen) {
    // Auto-focus input
    setTimeout(() => document.getElementById('aiChatInput').focus(), 250);
    // Set dynamic greeting based on analysis state
    updateAiGreeting();
  }
}

function updateAiGreeting() {
  const messages = document.getElementById('aiChatMessages');
  // Only update if only the default greeting exists
  if (messages.children.length !== 1) return;
  const greeting = currentData
    ? `Sequence analyzed. GC content: ${currentData.gc}%. Ask me anything about your DNA.`
    : 'Paste a DNA sequence and run an analysis — I\'ll help you interpret the results.';
  messages.children[0].querySelector('.ai-msg-bubble').textContent = greeting;
}

function addMessage(text, role) {
  const messages = document.getElementById('aiChatMessages');
  const div = document.createElement('div');
  div.className = `ai-msg ${role}`;
  div.innerHTML = `<div class="ai-msg-bubble">${text}</div>`;
  messages.appendChild(div);
  messages.scrollTop = messages.scrollHeight;
  return div;
}

async function sendAiMessage() {
  const input = document.getElementById('aiChatInput');
  const userMessage = input.value.trim();
  if (!userMessage) return;
  input.value = '';
  input.disabled = true;

  addMessage(userMessage, 'user');
  const typingDiv = addMessage('...', 'ai ai-typing');
  typingDiv.id = 'aiTyping';

  try {
    // Step 1: smartBioAI — instant local logic, no network call
    const logicReply = smartBioAI(userMessage, currentSeq, currentData);

    if (logicReply !== '__NO_ANSWER__') {
      console.log('LOGIC USED');
      document.getElementById('aiTyping')?.remove();
      const msgDiv = addMessage('', 'ai');
      await typewriterEffect(msgDiv.querySelector('.ai-msg-bubble'), logicReply, 18);
    } else {
      // Step 2: Groq API (free tier)
      console.log('GROQ USED');
      const groqReply = await askGroq(userMessage);
      document.getElementById('aiTyping')?.remove();
      const msgDiv = addMessage('', 'ai');
      await typewriterEffect(msgDiv.querySelector('.ai-msg-bubble'), groqReply, 12);
    }
  } catch (err) {
    document.getElementById('aiTyping')?.remove();
    addMessage('Error: ' + err.message, 'ai');
  } finally {
    input.disabled = false;
    input.focus();
  }
}

function typewriterEffect(element, text, speed = 15) {
  return new Promise(resolve => {
    const messages = document.getElementById('aiChatMessages');
    let i = 0;
    element.textContent = '';
    const interval = setInterval(() => {
      element.textContent += text[i];
      i++;
      messages.scrollTop = messages.scrollHeight;
      if (i >= text.length) {
        clearInterval(interval);
        resolve();
      }
    }, speed);
  });
}

// ── Ollama (DISABLED — runs locally only, requires: ollama serve && ollama pull llama3.2) ──
// async function askOllama(question) {
//   try {
//     const res = await fetch('http://localhost:11434/api/generate', {
//       method: 'POST',
//       headers: {'Content-Type': 'application/json'},
//       body: JSON.stringify({model: 'llama3.2', prompt: question, stream: true})
//     });
//     // ... streaming logic
//   } catch (err) {
//     return "Ollama is not running or model unavailable. Start it with: ollama serve";
//   }
// }

// ── Groq API (ACTIVE — free tier) ─────────────────────────────────────────
async function askGroq(question) {
  try {
    const res = await fetch('/ask_groq', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({
        question: question,
        sequence: currentSeq || ''
      })
    });
    const data = await res.json();
    if (data.answer) return data.answer;
    if (data.error) return `Error: ${data.error}`;
    return 'No response from Groq.';
  } catch (err) {
    return `Network error: ${err.message}. Make sure server is running at http://127.0.0.1:5001`;
  }
}

// Keep askAI as alias used by handleAskAi (AI Insights tab)
async function askAI(question) {
  const logicReply = smartBioAI(question, currentSeq, currentData);
  if (logicReply !== '__NO_ANSWER__') return logicReply;
  return await askGroq(question);
}

function smartBioAI(question, sequence, analysis) {
  const q = question.toLowerCase();
  const seq = sequence || currentSeq || '';
  const d = analysis || currentData;

  // GREETING
  if (q === 'hi' || q === 'hello' || q.includes('hey')) {
    return `Hi, I'm Armaani — your AI bioinformatics assistant.\n\nI can help you with:\n• GC content analysis\n• Nucleotide positions (e.g., 23rd base)\n• Base counts (A, T, G, C)\n• Restriction enzyme sites\n• ORF detection\n• Mutations and sequence insights\n\nAsk me anything about your DNA.`;
  }

  // POSITION — handles "position 5", "5th base", "23rd base", "what is at 10"
  const ordinalMatch = q.match(/(\d+)(?:st|nd|rd|th)?\s*(?:base|nucleotide|position)?/);
  if (q.includes('position') || q.includes('base') || q.includes('nucleotide')) {
    const match = q.match(/\d+/);
    if (match) {
      const pos = parseInt(match[0]) - 1;
      if (seq && pos >= 0 && pos < seq.length) {
        return `Base at position ${pos + 1} is "${seq[pos]}"`;
      }
      if (seq) return `Invalid position. Sequence is ${seq.length} bp (positions 1–${seq.length}).`;
    }
  }

  // BASE COUNT — "count A", "how many T", "number of G"
  if (q.includes('how many') || q.includes('count') || q.includes('number of')) {
    const base = q.match(/\b([atgc])\b/i);
    if (base) {
      const b = base[1].toUpperCase();
      const count = seq.split('').filter(x => x === b).length;
      return `There are ${count} "${b}" bases in your sequence.`;
    }
  }

  // GC CONTENT
  if (q.includes('gc') || q.includes('gc content')) {
    const gc = d ? parseFloat(d.gc) : null;
    if (gc === null) return 'Run an analysis first to get GC content.';
    if (gc < 40) return `GC content is ${gc}% — Low (AT-rich region).`;
    if (gc > 60) return `GC content is ${gc}% — High (stable DNA).`;
    return `GC content is ${gc}% — Moderate.`;
  }

  // DNA DEFINITION
  if (q.includes('what is dna') || q === 'dna') return 'DNA (Deoxyribonucleic Acid) is the molecule that carries the genetic instructions for the development, functioning, growth and reproduction of all known living organisms.';

  // BASE DEFINITIONS
  if (q.includes('what is a') || q === 'a') return 'A = Adenine, one of the four DNA bases. It pairs with Thymine (T).';
  if (q.includes('what is t') || q === 't') return 'T = Thymine, pairs with Adenine (A).';
  if (q.includes('what is g') || q === 'g') return 'G = Guanine, pairs with Cytosine (C). G-C bonds have 3 hydrogen bonds.';
  if (q.includes('what is c') || q === 'c') return 'C = Cytosine, pairs with Guanine (G).';

  // LENGTH
  if (q.includes('length') || q.includes('how long') || q.includes('how many bp')) {
    return `Sequence length is ${seq.length} base pairs.`;
  }

  // ENZYME
  if (q.includes('ecori'))   return 'EcoRI recognizes GAATTC and cuts between G and AATTC.';
  if (q.includes('bamhi'))   return 'BamHI recognizes GGATCC and cuts between G and GATCC.';
  if (q.includes('hindiii')) return 'HindIII recognizes AAGCTT and cuts between A and AGCTT.';
  if (q.includes('noti'))    return 'NotI recognizes GCGGCCGC — an 8-cutter, rare in most genomes.';
  if (q.includes('restriction') || q.includes('enzyme')) {
    const found = d ? Object.keys(d.restriction_sites || {}) : [];
    return found.length > 0
      ? `Restriction sites found: ${found.join(', ')}.`
      : 'No common restriction sites found. Use Lab Tools to run a digest.';
  }

  // ORF
  if ((q.includes('orf') || q.includes('open reading')) && !q.includes('what do you mean') && !q.includes('explain') && !q.includes('what is') && !q.includes('define')) {
    if (!d) return 'Run an analysis to detect ORFs.';
    const total = d.orfs_6frame ? d.orfs_6frame.reduce((s, f) => s + f.orfs.length, 0) : 0;
    return total > 0
      ? `${total} ORF(s) detected. ORFs start with ATG and end at a stop codon (TAA, TAG, TGA). Check the 6-Frame ORFs tab.`
      : 'No ORFs detected. The sequence may be non-coding or too short.';
  }

  // MUTATION
  if ((q.includes('mutation') || q.includes('snp')) && !q.includes('explain') && !q.includes('what do you mean') && !q.includes('what is')) {
    return 'Mutations can be:\n• Silent — no amino acid change\n• Missense — different amino acid\n• Nonsense — premature stop codon\n\nUse the Mutator tab to simulate one.';
  }

  // MELTING TEMP
  if ((q.includes('melt') || q.includes('tm') || q.includes('temperature')) && !q.includes('explain') && !q.includes('what is') && !q.includes('what do you mean')) {
    const tm = d ? d.tm : null;
    if (!tm) return 'Run an analysis to get the melting temperature.';
    return `Melting temperature is ${tm}°C. Recommended PCR annealing: ~${(parseFloat(tm) - 5).toFixed(1)}°C.`;
  }

  // PROTEIN
  if ((q.includes('protein') || q.includes('amino')) && !q.includes('explain') && !q.includes('what is') && !q.includes('what do you mean') && !q.includes('define')) {
    const pp = d ? d.protein_props : null;
    if (!pp) return 'No protein translated. The sequence may lack a start codon (ATG).';
    return `Protein: ${pp.length} amino acids, MW ${pp.mw.toLocaleString()} Da, pI ${pp.pi}.`;
  }

  // No match — signal Ollama fallback
  return '__NO_ANSWER__';
}

// Keep localAI as alias for backward compatibility
function localAI(question) {
  return smartBioAI(question, currentSeq, currentData);
}

// On load, always show input tab and hide analysis nav sections
document.addEventListener('DOMContentLoaded', () => {
  showTabByName('input');
  document.getElementById('analysisNav').style.display = 'none';
  document.getElementById('labNav').style.display = 'none';
  document.getElementById('insightsNav').style.display = 'none';
});
