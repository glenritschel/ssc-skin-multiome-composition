# One cell. Answers: WHICH filenames produced each hint, for the two series where
# the hint list fired without an obvious cause. Printing a boolean was not enough.
import re, urllib.request

HINTS = ['atac', 'multiome', 'peak', 'fragment', 'barcodes', 'matrix', 'features']

def names(gse):
    stem = gse[:-3] + 'nnn'
    base = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{gse}/suppl/'
    out = []
    try:
        html = urllib.request.urlopen(base, timeout=90).read().decode('utf8', 'ignore')
        out += [f for f in sorted(set(re.findall(r'href="([^"?][^"]*)"', html)))
                if not f.startswith('/') and not f.startswith('http')]
    except Exception as e:
        out.append(f'<<dir failed: {e}>>')
    try:
        fl = urllib.request.urlopen(base + 'filelist.txt', timeout=90).read().decode('utf8', 'ignore')
        for line in fl.splitlines()[1:]:
            c = line.split('\t')
            if len(c) > 1 and c[1].strip():
                out.append(c[1].strip())
    except Exception as e:
        out.append(f'<<filelist failed: {e}>>')
    return out

for gse in ['GSE195452']:
    fs = names(gse)
    print(gse, '--', len(fs), 'names total')
    for h in HINTS:
        m = [f for f in fs if h in f.lower()]
        if m:
            print(f'  hint {h!r}: {len(m)} file(s)')
            for f in m[:15]:
                print('     ', f)
            if len(m) > 15:
                print(f'      ... and {len(m)-15} more')
    print()
