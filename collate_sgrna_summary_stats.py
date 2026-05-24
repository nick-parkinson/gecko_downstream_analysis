# -*- coding: utf-8 -*-
"""
Collating CRISPR screen quality control statistics,
outputted by count_gecko_sgrnas.py script.

@author: nparkins
"""
import os
import sys
import pandas as pd

def collate_qc(folder, outname):
    stats = {}
    
    for file in os.listdir(folder):
        if file.startswith("pool"):
            poolname = "_".join(file.split("_")[1:-1])
            stats[poolname] = {}
            with open(os.path.join(folder,file)) as f:
                for line in f:
                    l = line.strip("\n").split(": ")
                    k = l[0]
                    dat = l[1]
                    if k == 'Total reads':
                        stats[poolname]['Total reads'] = float(dat)
                    elif 'rong length' in k:
                        v = int(float(dat.split(" ")[0]))
                        vp = round(float(dat.split("(")[1].strip("%)")),1)
                        stats[poolname]['Wrong length spacer (poss empty vector)'] = v
                        stats[poolname]['Wrong length %'] = vp
                    elif 'nmatched' in k:
                        v = int(float(dat.split(" ")[0]))
                        vp = round(float(dat.split("(")[1].strip("%)")),1)
                        stats[poolname]['Non-matching 20mer (poss library contamination)'] = v
                        stats[poolname]['Non-matching %'] = vp
                    elif 'perfect' in k:
                        v = int(float(dat.split(",")[0]))
                        vp = round(float(dat.split("(")[1].split("%")[0]),1)
                        stats[poolname]['Total perfectly mapped reads'] = v
                        stats[poolname]['% matched reads'] = vp
                    elif 'guides found' in k:
                        found = int(float(dat.split("/")[0]))
                        expected = int(dat.split("/")[1].split(" ")[0])
                        vp = round(float(dat.split("(")[1].strip("%)")),1)
                        stats[poolname]['Guides found'] = found
                        stats[poolname]['Expected library size'] = expected
                        stats[poolname]['% guides found'] = vp
                        stats[poolname]['Reads per sgRNA'] = stats[poolname]['Total perfectly mapped reads'] / expected
                    elif 'Skew' in k:
                        if dat.startswith('Not'):
                            stats[poolname]['Skew ratio'] = 'N/A'
                        else:
                            stats[poolname]['Skew ratio'] = float(dat)
                        
    summary = pd.DataFrame.from_dict(stats, orient='index')
    print(summary.head)
    try:
        summary.to_csv(os.path.join(folder,outname),
                       columns = ['Total reads',
                                   'Total perfectly mapped reads',
                                    '% matched reads',
                                     'Wrong length spacer (poss empty vector)',
                                     'Wrong length %',
                                     'Non-matching 20mer (poss library contamination)',
                                     'Non-matching %',
                                     'Guides found',
                                     'Expected library size',
                                      '% guides found',
                                      'Reads per sgRNA',
                                      'Skew ratio'
                                   ])
    except:
        summary.to_csv(os.path.join(folder,outname),
                       columns = ['Total reads',
                                   'Total perfectly mapped reads',
                                    '% matched reads',
                                     'Wrong length spacer (poss empty vector)',
                                     'Wrong length %',
                                     'Guides found',
                                     'Expected library size',
                                      '% guides found',
                                      'Reads per sgRNA',
                                      'Skew ratio'
                                   ])
        
if __name__ == "__main__":
    args = sys.argv[1:]
    collate_qc(args[0],args[1])
 