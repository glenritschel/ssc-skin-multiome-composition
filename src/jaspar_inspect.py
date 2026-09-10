# One cell. Ten seconds. Resolves what the pipeline actually scanned.
#
# The run reported: ['SMAD3','RUNX2','CREB1','EGR1','ELK4','HNF4A','HNF4A','POU5F1','Runx1','RUNX3']
# Two problems visible in that list alone:
#   HNF4A appears TWICE  -> two JASPAR matrices share the name; a dict keyed by name
#                           silently keeps one and drops the other.
#   'Runx1' is LOWERCASE -> JASPAR's legacy naming. Whether it is the human or the mouse
#                           matrix cannot be told from the name. Print the species.
from pyjaspar import jaspardb
import pandas as pd
motifs = jaspardb(release='JASPAR2024').fetch_motifs(collection='CORE',
                                                     tax_group=['vertebrates'])
WANT = ['RUNX1','RUNX2','RUNX3','SMAD3','EGR1','CREB1','ELK4','HNF4A','POU5F1']
rows = []
for m in motifs:
    if m.name.upper() in WANT:
        sp = getattr(m, 'species', None)
        rows.append(dict(name=m.name, matrix_id=m.matrix_id,
                         species=str(sp), length=m.length))
df = pd.DataFrame(rows).sort_values(['name','matrix_id'])
print(df.to_string(index=False))
print()
dup = df.name.str.upper().value_counts()
print('names with more than one matrix:', dict(dup[dup > 1]))
print()
print('9606 = Homo sapiens, 10090 = Mus musculus')
print('Any RUNX matrix that is not 9606 must be dropped or the paralog question is')
print('being answered with a mouse position weight matrix.')
