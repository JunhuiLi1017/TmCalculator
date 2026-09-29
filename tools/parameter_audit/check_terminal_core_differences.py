"""Reconstruct WO2001094611A2 Tables 2-3 increments from experimental Table 1.

Python standard library only. Values transcribed from printed pp47-54.
This checks the experimental-core baseline, not an unidentified fitting program.
All energies are kcal/mol; entropy columns are deliberately not reconstructed:
the separately reported experimental H/S/G columns need not be algebraically
identical after estimation/rounding (and some entries have larger discrepancies).
"""
import csv
from pathlib import Path

# Key, Table 1 core, full-duplex H, full-duplex G37, Tables 2-3 H, Tables 2-3 G37.
DATA = [
    ('AA/TA', 'TGAGCTCA', -56.7, -9.07, -3.1, -.67),
    ('TA/AA', 'AGAGCTCT', -55.1, -8.91, -2.5, -.58),
    ('CA/GA', 'GTAGCTAC', -60.3, -9.03, -4.3, -1.01),
    ('GA/CA', 'CGATATCG', -67.8, -8.87, -8.0, -.99),
    ('AC/TC', 'TGAGCTCA', -50.6, -8.15, -.1, -.21),
    ('TC/AC', 'AGAGCTCT', -51.4, -8.34, -.7, -.29),
    ('CC/GC', 'GTAGCTAC', -55.8, -8.04, -2.1, -.52),
    ('GC/CC', 'CGATATCG', -59.8, -8.12, -3.9, -.62),
    ('AG/TG', 'TGAGCTCA', -52.6, -8.57, -1.1, -.42),
    ('TG/AG', 'AGAGCTCT', -52.3, -8.34, -1.1, -.29),
    ('CG/GG', 'GTAGCTAC', -59.2, -8.68, -3.8, -.83),
    ('GG/CG', 'CGATATCG', -65.7, -8.81, -.7, -.96),
    ('AT/TT', 'TGAGCTCA', -55.4, -8.62, -2.4, -.45),
    ('TT/AT', 'AGAGCTCT', -56.5, -8.72, -3.2, -.48),
    ('CT/GT', 'GTAGCTAC', -63.8, -8.75, -6.1, -.87),
    ('GT/CT', 'CGATATCG', -66.8, -8.60, -7.4, -.86),
]
# Mixed-base cases: Table 1 pp48-51 -> Table 3 pp53-54.
DATA += [
    ('AA/TC', 'TGAGCTCA', -53.6, -8.42, -1.6, -.35),
    ('AC/TA', 'TGAGCTCA', -54.0, -8.92, -1.8, -.59),
    ('CA/GC', 'GTAGCTAC', -56.8, -8.53, -2.6, -.76),
    ('CC/GA', 'GTAGCTAC', -57.1, -8.71, -2.7, -.85),
    ('GA/CC', 'CGATATCG', -61.8, -8.30, -5.0, -.71),
    ('GC/CA', 'CGATATCG', -58.3, -8.91, -3.2, -1.01),
    ('TA/AC', 'AGAGCTCT', -54.6, -8.66, -2.3, -.45),
    ('TC/AA', 'AGAGCTCT', -55.5, -8.85, -2.7, -.55),
    ('AC/TT', 'TGAGCTCA', -52.2, -8.38, -.9, -.33),
    ('AT/TC', 'TGAGCTCA', -55.1, -8.42, -2.3, -.35),
    ('CC/GT', 'GTAGCTAC', -58.0, -8.39, -3.2, -.69),
    ('CT/GC', 'GTAGCTAC', -59.4, -8.21, -3.9, -.60),
    ('GC/CT', 'CGATATCG', -61.7, -8.33, -4.9, -.72),
    ('GT/CC', 'CGATATCG', -57.9, -8.11, -3.0, -.61),
    ('TC/AT', 'AGAGCTCT', -55.0, -8.80, -2.5, -.52),
    ('TT/AC', 'AGAGCTCT', -51.5, -8.44, -.7, -.34),
    ('AA/TG', 'TGAGCTCA', -54.2, -8.77, -1.9, -.52),
    ('AG/TA', 'TGAGCTCA', -55.4, -9.03, -2.5, -.65),
    ('CA/GG', 'GTAGCTAC', -59.4, -8.76, -3.9, -.88),
    ('CG/GA', 'GTAGCTAC', -63.7, -9.46, -6.0, -1.23),
    ('GA/CG', 'CGATATCG', -60.4, -8.50, -4.3, -.80),
    ('GG/CA', 'CGATATCG', -61.1, -9.04, -4.6, -1.08),
    ('TA/AG', 'AGAGCTCT', -54.0, -8.82, -2.0, -.53),
    ('TG/AA', 'AGAGCTCT', -54.8, -8.90, -2.4, -.57),
    ('AG/TT', 'TGAGCTCA', -56.8, -8.64, -3.2, -.45),
    ('AT/TG', 'TGAGCTCA', -57.4, -8.80, -3.5, -.54),
    ('CG/GT', 'GTAGCTAC', -59.2, -8.93, -3.8, -.96),
    ('CT/GG', 'GTAGCTAC', -64.8, -8.63, -6.6, -.81),
    ('GG/CT', 'CGATATCG', -63.3, -8.42, -5.7, -.76),
    ('GT/CG', 'CGATATCG', -63.7, -8.73, -5.9, -.92),
    ('TG/AT', 'AGAGCTCT', -57.8, -8.95, -3.9, -.59),
    ('TT/AG', 'AGAGCTCT', -57.3, -8.95, -3.6, -.59),
]

CORES = {
    'TGAGCTCA': (-50.5, -7.73),
    'AGAGCTCT': (-50.0, -7.76),
    'GTAGCTAC': (-51.6, -7.01),
    'CGATATCG': (-51.9, -6.89),
}

def main():
    rows = []
    for key, core, full_h, full_g, published_h, published_g in DATA:
        core_h, core_g = CORES[core]
        h, g = (full_h-core_h)/2, (full_g-core_g)/2
        rows.append(dict(key=key, core=core, full_H=full_h, core_H=core_h,
                         half_difference_H=round(h, 6), published_H=published_h,
                         full_G37=full_g, core_G37=core_g,
                         half_difference_G37=round(g, 6), published_G37=published_g,
                         H_matches_printed_precision=abs(h-published_h) <= .051,
                         G_matches_printed_precision=abs(g-published_g) <= .0051))
    assert len(rows) == 48
    assert all(row['G_matches_printed_precision'] for row in rows)
    assert [row['key'] for row in rows if not row['H_matches_printed_precision']] == ['GG/CG']
    path = Path(__file__).resolve().parents[2] / 'inst/extdata/terminal_core_reconstruction.tsv'
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    print('Tables 2-3: 48/48 G37 and 47/48 H match experimental-core half differences.')
    print('GG/CG: H half difference = -6.9, printed H = -0.7; no parameter correction made.')

if __name__ == '__main__':
    main()
