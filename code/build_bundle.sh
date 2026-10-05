#!/bin/bash
# Assemble the single deliverable zip: paper (pdf, md, tex), response to the referee, figures, code, data and logs.
# Usage: code/build_bundle.sh [output_dir] [referee_audit_dir]   (defaults: the parent of the project directory)
set -euo pipefail
G="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$(dirname "$G")}"; AUDIT="${2:-}"
NAME=Evidence_Press_IHC_heat_flux_bounds; STAGE="$(mktemp -d)"; B="$STAGE/$NAME"
mkdir -p "$B"/{paper,figs,code/cert,data/dns,data/steady,logs,problem}
cp "$G"/README.md "$G"/DEPENDENCIES.md "$G"/requirements.txt "$G"/reproduce.sh "$B"/
cp "$G"/paper/{main.pdf,paper.md,main.tex,abstract.tex,sec_*.tex,tab_*.tex,refs_apa.tex,refs_apa.md,references.bib,references.json,apa.csl,crossref.json,facts_*.json,novelty_report.md,response_to_referee.md,response_to_referee.pdf} "$B"/paper/
cp "$G"/figs/* "$B"/figs/
cp "$G"/code/*.py "$G"/code/*.sh "$G"/code/*.json "$G"/code/*.pkl "$G"/code/*.log "$B"/code/
cp "$G"/code/cert/*.py "$B"/code/cert/; mkdir -p "$B"/code/superseded && cp "$G"/code/superseded/* "$B"/code/superseded/
cp "$G"/data/dns/* "$B"/data/dns/
cp "$G"/data/steady/* "$B"/data/steady/
cp "$G"/logs/* "$B"/logs/
cp "$G"/Open_problems_geology_frontier_AI.md "$B"/problem/
if [[ -n "$AUDIT" && -d "$AUDIT" ]]; then mkdir -p "$B"/referee_audit && cp -r "$AUDIT"/* "$B"/referee_audit/; fi
# The field snapshots of the simulations (data/dns/*_snap.npz) go into a companion zip, so that each archive stays
# below 30 MiB. Nothing in verify_numbers.py reads them; only dns_resolution.py does.
S2="$(mktemp -d)"; mkdir -p "$S2/$NAME/data/dns"; mv "$B"/data/dns/*_snap.npz "$S2/$NAME/data/dns/"
( cd "$B" && find . -type f ! -name MANIFEST.sha256 -print0 | sort -z | xargs -0 sha256sum > MANIFEST.sha256 )
( cd "$S2/$NAME" && find . -type f ! -name 'MANIFEST*' -print0 | sort -z | xargs -0 sha256sum > MANIFEST_snapshots.sha256 )
rm -f "$OUT/$NAME.zip" "$OUT/${NAME}_snapshots.zip"
( cd "$STAGE" && zip -9 -qr "$OUT/$NAME.zip" "$NAME" )
( cd "$S2" && zip -9 -qr "$OUT/${NAME}_snapshots.zip" "$NAME" )
rm -rf "$STAGE" "$S2"
for z in "$OUT/$NAME.zip" "$OUT/${NAME}_snapshots.zip"; do ls -la "$z"; unzip -l "$z" | tail -1; done
