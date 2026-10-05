"""Build the APA 7 reference list for the LaTeX paper (natbib-compatible thebibliography) and a CSL-JSON file for
the Markdown version. Formatting is delegated to pandoc's citeproc with the official APA 7 CSL style."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import tempfile, json, subprocess, re, os
TMPF = os.path.join(tempfile.mkdtemp(), 'nocite.md')
P = (ROOT + '/paper/')
csl = json.loads(subprocess.run(['pandoc', P+'references.bib', '-t', 'csljson'], capture_output=True, text=True, check=True).stdout)
for e in csl:   # citeproc title-cases container titles; keep Italian particles lower case
    if e['id'] == 'Talenti1976': e['container-title'] = 'Annali <span class="nocase">della</span> Scuola Normale Superiore <span class="nocase">di</span> Pisa, Classe <span class="nocase">di</span> Scienze (4)'
for e in csl:
    for a in e.get('author', []):   # APA keeps Dutch particles with the surname: "van der Poel, E. P."
        if 'dropping-particle' in a: a['non-dropping-particle'] = a.pop('dropping-particle')
    if '-' in str(e.get('issue', '')): e['issue'] = str(e['issue']).replace('-', '–')
json.dump(csl, open(P+'references.json', 'w'), indent=1, ensure_ascii=False)
open(TMPF, 'w').write("---\nnocite: '@*'\n---\n")
md = subprocess.run(['pandoc', TMPF, '--citeproc', '--bibliography', P+'references.json', '--csl', P+'apa.csl',
                     '-t', 'markdown-citations-smart', '--wrap=none'], capture_output=True, text=True, check=True).stdout
md = re.sub(r'\[([^\]]+)\]\{\.nocase\}', r'\1', md)
entries = re.findall(r'::: \{#ref-(\S+) \.csl-entry\}\n(.*?)\n:::', md, flags=re.S)
by_id = {e['id']: e for e in csl}
def label(e):
    a = [x.get('family', x.get('literal', '')) for x in e['author']]
    a = [(x.get('non-dropping-particle', '') + ' ' + f).strip().replace(' ', '~') if x.get('non-dropping-particle') else f
         for x, f in zip(e['author'], a)]
    y = str(e['issued']['date-parts'][0][0])
    short = a[0] if len(a) == 1 else (f'{a[0]} and {a[1]}' if len(a) == 2 else f'{a[0]} et~al.')
    full = ', '.join(a[:-1]) + ', and ' + a[-1] if len(a) > 2 else short
    return f'{short}({y}){full}'
out = [r'\begin{thebibliography}{%d}' % len(entries), r'\providecommand{\natexlab}[1]{#1}', r'\setlength{\itemsep}{2pt}']
for key, text in entries:
    tex = subprocess.run(['pandoc', '-f', 'markdown-smart', '-t', 'latex', '--wrap=none'], input=text, capture_output=True, text=True, check=True).stdout.strip()
    out.append(f'\\bibitem[{label(by_id[key])}]{{{key}}}\n{tex}')
out.append(r'\end{thebibliography}')
open(P+'refs_apa.tex', 'w').write('\n\n'.join(out) + '\n')
open(P+'refs_apa.md', 'w').write('\n\n'.join(t for _, t in entries) + '\n')
print(len(entries), 'entries')
