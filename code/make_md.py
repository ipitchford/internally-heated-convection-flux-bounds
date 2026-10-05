"""Markdown version of the paper. Flattens the LaTeX sources, resolves theorem, equation, figure, table and
section numbers from main.aux (so that the Markdown matches the PDF), points figures at the PNG renderings and
converts with pandoc, rendering APA 7 citations and references via citeproc and the official APA CSL style."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import tempfile
TMPF = os.path.join(tempfile.mkdtemp(), 'main_flat.tex')
import re, subprocess
P = (ROOT + '/paper/')
aux = dict(re.findall(r'\\newlabel\{([^}]+)\}\{\{([^}]*)\}', open(P+'main.aux').read()))
def flatten(fn):
    s = open(P+fn).read()
    return re.sub(r'\\input\{([^}]+)\}', lambda m: flatten(m.group(1) if m.group(1).endswith('.tex') else m.group(1)+'.tex'), s)
s = flatten('main.tex')
s = re.sub(r'(?m)^%.*$', '', s)
s = s.replace(r'\input{refs_apa.tex}', '')
# theorem-like environments -> numbered bold run-in headings that match the PDF
names = dict(theorem='Theorem', proposition='Proposition', lemma='Lemma', corollary='Corollary', remark='Remark')
sec, cnt, out, pos = 0, 0, [], 0
pat = re.compile(r'\\section\{|\\appendix|\\begin\{(theorem|proposition|lemma|corollary|remark)\}(\[([^\]]*)\])?(\s*\\label\{([^}]*)\})?')
for m in pat.finditer(s):
    out.append(s[pos:m.start()]); pos = m.end(); t = m.group(0)
    if t.startswith('\\section{'): sec += 1; cnt = 0; out.append(t); continue
    if t == '\\appendix': out.append(t); continue
    cnt += 1; num = f'{sec}.{cnt}'
    if m.group(5): assert aux[m.group(5)] == num, (m.group(5), aux[m.group(5)], num)
    head = f'\\textbf{{{names[m.group(1)]} {num}}}' + (f' ({m.group(3)})' if m.group(3) else '') + '. '
    out.append('\n\n' + head)
out.append(s[pos:]); s = ''.join(out)
s = re.sub(r'\\end\{(theorem|proposition|lemma|corollary|remark)\}', '\n\n', s)
s = re.sub(r'\\begin\{proof\}\[([^\]]*)\]', lambda m: '\n\n\\emph{' + m.group(1) + '.} ', s)
s = s.replace(r'\begin{proof}', '\n\n\\emph{Proof.} ').replace(r'\end{proof}', ' $\\square$\n\n')
# cross-references and equation tags from the aux file
s = re.sub(r'\\eqref\{([^}]+)\}', lambda m: '(' + aux[m.group(1)] + ')', s)
s = re.sub(r'\\ref\{([^}]+)\}', lambda m: aux[m.group(1)], s)
s = re.sub(r'\\label\{(eq:[^}]+)\}', lambda m: '\\tag{' + aux[m.group(1)] + '}', s)
def number_captions(s):
    """Prefix each float caption with 'Figure N.' or 'Table N.' using the label that follows it."""
    out, i = [], 0
    while True:
        j = s.find('\\caption{', i)
        if j < 0: out.append(s[i:]); break
        k, depth = j + len('\\caption{'), 1
        while depth:
            depth += {'{': 1, '}': -1}.get(s[k], 0); k += 1
        body = s[j+len('\\caption{'):k-1]
        lab = re.search(r'\\label\{((fig|tab):[^}]+)\}', body + s[k:k+80])
        pre = ('Figure ' if lab.group(2) == 'fig' else 'Table ') + aux[lab.group(1)] + '. ' if lab else ''
        out.append(s[i:j] + '\\caption{' + pre + body + '}'); i = k
    return ''.join(out)
s = number_captions(s)
s = re.sub(r'\\includegraphics(\[[^\]]*\])?\{([^}]+)\.pdf\}', r'\\includegraphics{../figs/\2.png}', s)
s = s.replace(r'\appendix', 'REFSPLACEHOLDER\n\n\\appendix')
open(TMPF, 'w').write(s)
r = subprocess.run(['pandoc', TMPF, '-f', 'latex', '-t', 'markdown-citations-smart-raw_attribute-raw_html-raw_tex-simple_tables-multiline_tables-grid_tables+pipe_tables+tex_math_dollars',
                    '--citeproc', '--bibliography', P+'references.json', '--csl', P+'apa.csl', '--wrap=none', '-s'],
                   capture_output=True, text=True)
print(r.returncode, r.stderr[-2000:])
md = r.stdout
# move the reference list to where the PDF has it (before the appendices) and tidy the metadata
idx = md.find('\n::: {#refs')
refs, body = (md[idx:].rstrip() + '\n', md[:idx]) if idx >= 0 else ('', md)
body = body.replace('REFSPLACEHOLDER', '# References {.unnumbered}\n' + refs)
body = re.sub(r'(?m)^(bibliography|csl|link-citations):.*\n', '', body)
body = re.sub(r'\[\]\{#[^}]*\}', '', body)
body = re.sub(r'\[([^\]]+)\]\{\.nocase\}', r'\1', body)
open(P+'paper.md', 'w').write(body)
print(len(body), 'chars;', body.count('csl-entry'), 'reference entries;', len(re.findall(r'\]\(\.\./figs/', body)), 'figures')
# section numbers matching the PDF (numbered sections; lettered appendices after the reference list)
lines, s1, s2, app = body.split('\n'), 0, 0, False
for i, l in enumerate(lines):
    if l.startswith('# ') and '.unnumbered' not in l:
        if app: s1 = 'A' if s1 is None else chr(ord(s1)+1)
        else: s1 += 1
        s2 = 0; lines[i] = f'# {s1} ' + l[2:]
    elif l.startswith('## ') and not app:
        s2 += 1; lines[i] = f'## {s1}.{s2} ' + l[3:]
    elif l.startswith('# References'):
        app = True; s1 = None
body = '\n'.join(lines)
open(P+'paper.md', 'w').write(body)
