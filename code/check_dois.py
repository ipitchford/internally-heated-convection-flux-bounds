"""Check every reference against Crossref. For each entry of paper/references.bib that has a DOI, fetch the
Crossref record of that DOI and compare its year, first-author surname and title with the BibTeX entry. Entries
without a DOI are listed with their URL. Writes paper/crossref.json; exits with status 1 on any mismatch."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import json, re, time, unicodedata, urllib.request, urllib.parse


def norm(s):
    """lower-case ASCII words of a title, with LaTeX accents, braces and HTML tags removed"""
    s = re.sub(r'<[^>]+>', ' ', s)
    s = re.sub(r"\\[`'^\"~=.uvHc]\s*", '', s)          # accent commands: {\'e} -> {e}
    s = s.replace('{', '').replace('}', '').replace('&amp;', '&')
    for dash in ('--', '\u2013', '\u2014', '\u2010', '-'): s = s.replace(dash, ' ')
    s = unicodedata.normalize('NFKD', s).encode('ascii', 'ignore').decode().lower()
    return re.sub(r'[^a-z0-9]+', ' ', s).split()


def get(url, tries=3):
    req = urllib.request.Request(url, headers={'User-Agent': 'evidence-press-reference-check/1.0'})
    for i in range(tries):
        try:
            return json.load(urllib.request.urlopen(req, timeout=30))
        except Exception:
            if i == tries-1: raise
            time.sleep(3*(i+1))


bib = open(ROOT + '/paper/references.bib').read()
entries = {}
for m in re.finditer(r'@(\w+)\{(\w+),(.*?)\}\s*(?=@|\Z)', bib, re.S):
    fields = dict((k.lower(), v.strip()) for k, v in re.findall(r'(\w+)\s*=\s*\{((?:[^{}]|\{[^{}]*\})*)\}', m.group(3)))
    entries[m.group(2)] = fields
out, bad = {}, []
for key, f in entries.items():
    rec = dict(bib_title=f.get('title'), bib_year=f.get('year'), bib_doi=f.get('doi'), bib_url=f.get('url'))
    if f.get('doi'):
        try:
            m = get('https://api.crossref.org/works/' + urllib.parse.quote(f['doi']))['message']
            rec.update(title=(m.get('title') or [''])[0], journal=(m.get('container-title') or [''])[0], vol=m.get('volume'),
                       issue=m.get('issue'), page=m.get('page') or m.get('article-number'),
                       year=(m.get('issued') or m.get('published') or {'date-parts': [[None]]})['date-parts'][0][0],
                       year_print=(m.get('published-print') or {'date-parts': [[None]]})['date-parts'][0][0],
                       authors=[a.get('family', '') + ', ' + a.get('given', '') for a in m.get('author', [])])
            bt, ct = set(norm(f.get('title', ''))), set(norm(rec['title']))
            first = norm(f.get('author', '').split(' and ')[0].split(',')[0])
            cfirst = norm(rec['authors'][0].split(',')[0]) if rec['authors'] else []
            rec['title_overlap'] = len(bt & ct)/max(1, len(bt))
            rec['match'] = bool(rec['title_overlap'] >= 0.8 and f.get('year') in (str(rec['year']), str(rec['year_print']))
                                and (not first or not cfirst or first[-1] == cfirst[-1]))
        except Exception as e:
            rec.update(error=str(e), match=False)
        time.sleep(0.3)
    else:
        rec['match'] = None
    out[key] = rec
    if rec['match'] is False: bad.append(key)
    print(f"{key:28s} {'OK ' if rec['match'] else ('-- ' if rec['match'] is None else 'BAD')} "
          f"{rec.get('year', rec['bib_year'])} | {(rec.get('title') or rec['bib_title'] or '')[:70]} | {rec.get('journal', '')} "
          f"{rec.get('vol') or ''} {rec.get('page') or ''} | {f.get('doi') or f.get('url') or 'no DOI/URL'}")
json.dump(out, open(ROOT + '/paper/crossref.json', 'w'), indent=1)
print(f"\n{sum(1 for r in out.values() if r['match'])} of {sum(1 for r in out.values() if r['match'] is not None)} DOIs match their "
      f"BibTeX entries; {sum(1 for r in out.values() if r['match'] is None)} entries have no DOI.")
if bad:
    print('MISMATCH:', ', '.join(bad)); sys.exit(1)
