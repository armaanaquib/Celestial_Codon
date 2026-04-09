# Celestial Codon — DNA Analysis Platform

A professional bioinformatics web application for analyzing DNA sequences, built with Flask and powered by the Groq AI API.

## Features

- **Sequence Analysis** — GC content, AT content, melting temperature, nucleotide frequency
- **6-Frame ORF Detection** — Canvas visualization of all 3 forward + 3 reverse reading frames
- **Protein Analysis** — Molecular weight, isoelectric point (pI), hydrophobicity profile
- **Codon Bias** — Compare codon usage against E. coli, Human, and Yeast reference tables
- **Restriction Enzymes** — 30+ enzymes with N-wildcard support
- **Motif Library** — TATA-box, Kozak, Shine-Dalgarno, NF-kB, SP1, AP1 and more
- **Mutation Simulator** — Point mutation with synonymous/missense/nonsense classification
- **Disease Lookup** — MyVariant.info + ClinVar + dbSNP integration
- **Lab Tools** — Primer design, restriction digest simulator, PCR product predictor
- **Secondary Structure** — CpG island detection + Chou-Fasman prediction
- **AI Assistant (Armaani)** — Powered by Groq API (llama-3.1-8b-instant), free tier
- **Export** — Excel (.xlsx with colored sheets), FASTA, JSON

## Setup

### 1. Install dependencies
```bash
pip install -r requirements.txt
```

### 2. Get a free Groq API key
Sign up at https://console.groq.com and create an API key.

### 3. Set environment variable
```bash
# Windows
setx GROQ_API_KEY "your_key_here"

# Linux / Mac
export GROQ_API_KEY="your_key_here"
```

### 4. Run
```bash
python app.py
```

Open http://127.0.0.1:5001 in your browser.

## Deployment (Render)

1. Push this repo to GitHub
2. Go to https://render.com → New Web Service → connect your repo
3. Set build command: `pip install -r requirements.txt`
4. Set start command: `gunicorn app:app`
5. Add environment variable: `GROQ_API_KEY` = your key
6. Deploy

## Project Structure

```
celestial-codon/
├── app.py              # Flask backend — all routes and analysis logic
├── static/
│   ├── app.js          # Frontend JavaScript — rendering, AI chat, exports
│   ├── style.css       # Dark theme stylesheet
│   ├── armaan2.png     # AI assistant avatar
│   └── style_test.txt  # (ignore)
├── templates/
│   └── index.html      # Single-page app HTML
├── requirements.txt
├── Procfile            # For Render/Heroku deployment
└── README.md
```

## Tech Stack

- **Backend**: Python, Flask, Matplotlib, NumPy, openpyxl
- **Frontend**: Vanilla JavaScript, CSS custom properties
- **AI**: Groq API (Meta Llama 3.1 8B Instant) — free tier
- **Charts**: Matplotlib rendered server-side, base64 encoded

## License

Open source — feel free to use, modify, and build on this project.
