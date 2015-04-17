"""
Read Granular Layer netowrk structure from files
generate by the Solinas2010 model.
"""
import numpy as np
# import ipdb

def read_lol_file(dirname, filename,separator = None):
    if filename is None:
        print 'The required filemane is missing.'
        return None
    try:
        f = open(dirname+filename)
    except IOError:
        print 'Cannot open', dirname+filename
        
    if separator is None:
        lol = [[int(t) for t in l[:-1].split()] for l in f.readlines() ]
    else:
        lol = [[int(t) for t in l[:-1].split(separator)] for l in f.readlines()]
    f.close()
    return lol

def read_matrix_file(dirname, filename,separator = None):
    if filename is None:
        print 'The required filemane is missing.'
        return None
    try:
        matrix = np.loadtxt(dirname+filename)
    except IOError:
        print 'Cannot open', dirname+filename
        
    return matrix

def dig_structure_dict(data_dir,d):
    for it_n,it in d.iteritems():
        if 'filename' in it.keys():
            print 'loading: ', it['filename'], it['Comment'], '...'
            if it_n == 'positions':
                it['data'] = read_matrix_file(data_dir,it['filename'])
            elif it['filename'] is None:
                it['data'] = None
            else:
                it['data'] = read_lol_file(data_dir,it['filename'])
        else:
            dig_structure_dict(data_dir,it)
    return None

def load_Solinasetal2010_structure():
    data_dir = 'SimData/'
    structure = {'Glomeruli':
                 {'positions':{'filename':'glom_coord.lst','Comment':'3D position of glomeruli: (x, y, z, idx)'},
                  'convergence':{'filename':None,'Comment':'Convergence of Golgi cell axons into glomeruli'},
                  'divergence_to_grc':{'filename':'div_glom_target_grc.lst','Comment':'Divergence of Glomruli to granule cells'},
                  'divergence_to_goc':{'filename':'div_glom_target_goc.lst','Comment':'Divergence of Glomruli to Golgi cells'}
              },
                 'Golgis':{'positions':{'filename':'goc_coord.lst','Comment':'3D position of Golgi cells: (x, y, z, idx)'},
                  'convergence':{'filename':'conv_goc_glom_sources.lst','Comment':'Convergence of glimeruli to a single Golgi cell'},
                  'divergence_to_glom':{'filename':'div_goc_targets.lst','Comment':'Divergence of Golgi cell axons to Glomeruli'}},
                 'Granules':{'positions':{'filename':'grc_coord.lst','Comment':'3D position of granule cells: (x, y, z, idx)'},
                  'convergence':{'filename':'conv_goc_glom_sources.lst','Comment':'Convergence of glomeruli to Golgi cells'},
                  'divergence_to_goc':{'filename':'div_grc_targets.lst','Comment':'Divergence of granule cell axons to Golgi cells'}}
             }

    print structure

    dig_structure_dict(data_dir,structure)

    return structure


if __name__ == "__main__":
    load_Solinasetal2010_structure()
