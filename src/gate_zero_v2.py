# GATE ZERO v2 --- run this as one Colab cell before anything else.
#
# Fixes a defect in v1: v1 listed the FTP suppl/ directory only, so a series
# whose per-sample files sit inside GSE*_RAW.tar looked empty. GEO publishes
# filelist.txt, which enumerates the tar contents. v1 never read it, so it
# could return a false NONE. This version reads it.

import re, urllib.request

CANDIDATES = {
    'GSE195452': 'Gur 2022 Cell, LGR5 fibroblast hub. Believed scRNA + CITE-seq.',
    'GSE312129': 'JCI Insight 2025 (PMID 41411065). PAIRED snRNA + snATAC multiome, '
                 '10 dcSSc vs 4 HC, lesional forearm.',
    'GSE312932': 'Same study, spatial arm.',
    'GSE99702' : 'Liu 2020 Nat Commun. BULK ATAC on sorted skin cell types.',
}
HINTS = ['atac', 'multiome', 'peak', 'fragment', 'tbi', 'barcodes', 'matrix', 'features']

def listing(gse):
    stem = gse[:-3] + 'nnn'
    base = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{gse}/suppl/'
    out = []
    try:
        html = urllib.request.urlopen(base, timeout=90).read().decode('utf8', 'ignore')
        out += [f for f in sorted(set(re.findall(r'href="([^"?][^"]*)"', html)))
                if not f.startswith('/') and not f.startswith('http')]
    except Exception as e:
        out.append(f'<<dir listing failed: {e}>>')
    # the part v1 missed
    try:
        fl = urllib.request.urlopen(base + 'filelist.txt', timeout=90).read().decode('utf8', 'ignore')
        for line in fl.splitlines()[1:]:
            cols = line.split('\t')
            if len(cols) > 1 and cols[1].strip():
                out.append('  [in RAW.tar] ' + cols[1].strip())
    except Exception as e:
        out.append(f'  <<filelist.txt unavailable: {e}>>')
    return out

verdict = {}
for gse, note in CANDIDATES.items():
    print('=' * 74)
    print(gse, '--', note)
    print('=' * 74)
    files = listing(gse)
    for f in files[:60]:
        print('  ', f)
    if len(files) > 60:
        print(f'   ... and {len(files) - 60} more')
    blob = ' '.join(files).lower()
    hit = [h for h in HINTS if h in blob]
    verdict[gse] = hit
    print('  hints:', hit or 'NONE')
    print()

print('=' * 74)
for gse, hit in verdict.items():
    print(f'{gse:<12} {"CHROMATIN-LIKE FILES PRESENT" if hit else "no chromatin-like files"}')
print('=' * 74)

if not verdict.get('GSE312129'):
    raise RuntimeError(
        'GATE ZERO v2 FAILED for GSE312129, which is the only candidate with PAIRED\n'
        'RNA and ATAC in the same nuclei. Without pairing the paralog attribution step\n'
        'is impossible and the study reduces to replicating an already published result.\n'
        'Do not proceed. Report what printed above.'
    )
print('\nGATE ZERO v2 PASSED for GSE312129. Proceed.')
