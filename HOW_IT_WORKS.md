# How Celestial Codon Works

This document explains the full architecture, data flow, and every major component of the project.

---

## Overview

Celestial Codon is a **single-page web application** built with:
- **Flask** (Python) as the backend server
- **Vanilla JavaScript** as the frontend (no React, no Vue — just plain JS)
- **Matplotlib** for generating charts server-side
- **Groq API** for the AI assistant (free tier)
- **openpyxl** for Excel export

When a user pastes a DNA sequence and clicks Analyze, the browser sends it to the Flask server, which runs all the biology calculations and returns a big JSON object. The JavaScript then renders everything into the UI.

---

## File Structure

```
app.py              — The entire backend (Flask routes + all biology logic)
templates/
  index.html        — The entire frontend HTML (one page, many tabs)
static/
  app.js            — All frontend JavaScript (rendering, AI chat, exports)
  style.css         — All CSS (dark theme, layout, components)
  armaan2.png       — AI assistant avatar image
requirements.txt    — Python dependencies
Procfile            — Tells Render/Heroku how to start the app
```

---

## Backend: app.py

### How a DNA analysis works (step by step)

1. User sends `POST /analyze` with JSON `{ sequence, organism, theme, motif }`
2. `validate(seq)` — strips spaces/newlines, checks only A/T/C/G, min 3 bases
3. All calculations run in sequence:
   - `gc_content()` — counts G+C divided by total length
   - `melting_temp()` — Wallace rule for short sequences (<14bp), Marmur formula for longer
   - `reverse_complement()` — translates ATCG→TAGC then reverses
   - `transcribe()` — replaces T with U to get mRNA
   - `translate_seq()` — reads codons from AUG start, stops at UAA/UAG/UGA
   - `find_all_orfs_6frame()` — scans all 3 forward + 3 reverse frames for ATG→stop regions
   - `protein_properties()` — calculates MW (sum of AA weights), pI (binary search on charge curve), hydrophobicity (Kyte-Doolittle scale)
   - `codon_bias_comparison()` — counts codons per 1000, compares to hardcoded E.coli/Human/Yeast tables
   - `restriction_sites()` — searches for 30 enzyme recognition sequences, supports N wildcard via regex
   - `palindromes()` — finds sequences that equal their own reverse complement (length 4-12)
   - `scan_motif_library()` — regex search for TATA-box, Kozak, Shine-Dalgarno, NF-kB, etc.
   - `design_primers()` — takes first/last 20bp, calculates Tm and GC%, checks for hairpin risk
   - `make_charts()` — generates Matplotlib figures, saves to BytesIO buffer, returns base64 string
   - `make_circular_plot()` — plasmid-style circular map with GC arc coloring
4. All results packed into one JSON object and returned to the browser

### Other routes

| Route | What it does |
|---|---|
| `POST /mutate` | Simulates a point mutation, classifies as synonymous/missense/nonsense |
| `POST /disease` | Queries MyVariant.info + NCBI ClinVar for disease associations |
| `POST /digest` | Simulates restriction enzyme digest, returns fragment sizes |
| `POST /pcr` | Checks if two primers amplify a product from the sequence |
| `POST /secondary` | CpG island detection + Chou-Fasman secondary structure prediction |
| `POST /eli5` | Generates plain-English explanations of the analysis results |
| `POST /quiz` | Generates 5 multiple-choice questions based on the sequence |
| `POST /ask_groq` | Sends a question to Groq API (llama-3.1-8b-instant), returns answer |
| `POST /export/xlsx` | Generates a colored Excel file with 6 sheets using openpyxl |
| `POST /export/fasta` | Returns sequence in FASTA format |
| `POST /export/json` | Returns full analysis as a JSON file |

### Biology reference data (hardcoded in app.py)

- `CODON_TABLE` — all 64 RNA codons mapped to amino acids
- `AA_MW` — molecular weights of all 20 amino acids (Daltons)
- `AA_PKA` — pKa values for ionizable side chains (used for pI calculation)
- `AA_HYDRO` — Kyte-Doolittle hydrophobicity scale
- `RESTRICTION_ENZYMES` — 30 enzyme recognition sequences
- `CODON_BIAS` — codon usage frequency tables for E.coli K12, Human, S.cerevisiae
- `MOTIF_LIBRARY` — known promoter/TF binding motifs with degenerate base support

---

## Frontend: app.js

### How the UI works

The page has many tabs but only one is visible at a time. After analysis:
1. `analyze()` sends the sequence to `/analyze`
2. On success, `renderAll(data)` is called
3. `renderAll()` calls each render function for every tab
4. The sidebar navigation sections (Analysis, Lab Tools, Insights) become visible
5. The Overview tab is shown automatically

### Render functions

| Function | What it renders |
|---|---|
| `renderOverview(d)` | Stats cards, color-coded sequence viewer, freq bars, restriction table |
| `renderORFs(d)` | Canvas drawing of 6 reading frames + ORF list below |
| `renderProtein(d)` | MW/pI/hydrophobicity stats + domain map (red=hydrophobic, blue=hydrophilic) |
| `renderCodons(d)` | Codon frequency table with diff bars (green=higher, red=lower than reference) |
| `renderMotifLib(d)` | Motif hits or improved empty state with suggestions |
| `renderCharts(d)` | Sets `<img src>` to base64 chart from server |
| `renderAIInsights(d)` | Rule-based insight cards (GC stability, ORF count, protein type, etc.) |

### AI Assistant flow

```
User types message
       ↓
smartBioAI(question, sequence, data)
  — Checks for known patterns: position, base count, GC, length, enzymes, ORFs, etc.
  — Returns instant answer OR "__NO_ANSWER__" if it doesn't know
       ↓
If __NO_ANSWER__:
  askGroq(question)
  — POST /ask_groq → Flask → Groq API → llama-3.1-8b-instant
  — Returns answer as text
       ↓
typewriterEffect(bubble, text, speed)
  — Reveals text character by character like ChatGPT
```

### Why smartBioAI exists

Groq has rate limits. For simple questions like "what is the GC content" or "how many A bases", there's no need to call an API — the answer is already in `currentData`. `smartBioAI` handles these instantly with zero latency and zero API usage.

---

## Styling: style.css

Uses CSS custom properties (variables) for the entire color system:

```css
--bg          #0b1220   Main background (deep navy)
--surface     #131f35   Card background
--border      #1e2d47   Card borders
--accent      #3b82f6   Blue — buttons, active states
--green       #22c55e   GC content
--orange      #f97316   AT content
--purple      #a855f7   Melting temperature
--red         #ef4444   Warnings, errors
```

Light theme is a full override via `[data-theme=light]` on the `<html>` element.

---

## AI Assistant: Armaani

- Name: Armaani
- Avatar: `static/armaan2.png`
- Model: `llama-3.1-8b-instant` via Groq API
- Training cutoff: early 2024
- Fallback: `smartBioAI` local logic (no API needed)
- The assistant is bioinformatics-focused but can answer general questions too

To change the model, edit this line in `app.py`:
```python
model="llama-3.1-8b-instant"
```
Other free Groq models: `llama-3.3-70b-versatile`, `mixtral-8x7b-32768`

---

## Excel Export

The `/export/xlsx` route generates a `.xlsx` file with 6 sheets:
1. **Summary** — basic stats, sequences, protein properties
2. **ORFs** — all open reading frames across 6 frames
3. **Restriction & Palindromes** — enzyme cut sites and palindromic sequences
4. **Codon Bias vs E.coli** — codon frequency comparison (green=higher, red=lower)
5. **Amino Acids** — amino acid composition of the translated protein
6. **Charts** — pie chart of nucleotide composition

Column widths and row heights are calculated dynamically based on content length.

---

## Deployment

The app is designed to run on any platform that supports Python:
- **Local**: `python app.py` → http://127.0.0.1:5001
- **Render** (recommended free hosting): uses `Procfile` with `gunicorn app:app`
- **Railway / Heroku**: same Procfile works

The only required environment variable is `GROQ_API_KEY`.
