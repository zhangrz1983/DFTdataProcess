

from pymongo import MongoClient
import pprint

db = MongoClient()
search_dir = 'stacking_fault_and_dislocation/coarse_relax_stacking_fault'
search_crit = {'dir_name': {'$regex': search_dir}}
docs = list(db.vasp.tasks.find(search_crit))

doc_noshift = []
ene = {}
for doc in docs:
    if doc['state'] != 'killed':
        #ene[doc['dir_name'].split('/')[-1]] = doc['output']['final_energy']
        dir = doc['dir_name'].split('/')[-1].split('_')
        if len(dir) == 5:
            mater = dir[0]
            slip_sys = dir[1] + '_' + dir[3]
            fault_shift = float(dir[-1])
            if mater not in ene:
                ene[mater] = {}
            if slip_sys not in ene[mater]:
                ene[mater][slip_sys] = []
            ene[mater][slip_sys].append([fault_shift, doc['output']['final_energy']])
        elif len(dir) == 3:
            doc_noshift.append(doc)
        else:
            pass

fault_shift = 0.0
for doc in doc_noshift:
    dir = doc['dir_name'].split('/')[-1].split('_')
    mater = dir[0]
    plane = dir[1]
    for slip_sys in ene[mater].keys():
        if plane == slip_sys.split('_')[0]:
            ene[mater][slip_sys].append([fault_shift, doc['output']['final_energy']])
            ene[mater][slip_sys].sort()
            for e in ene[mater][slip_sys]:
                e[1] = round(e[1] - doc['output']['final_energy'], 3)

with open('stacking_fault_energy', 'w') as f:
    pprint.pprint(ene, f)

