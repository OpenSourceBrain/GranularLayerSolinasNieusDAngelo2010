"""
PYTHON script to run the network structure generation 
used in the Solinas et al. 2010
and load it into a PYTHON dictionary
"""
import load_structure as loader
import os
import ipdb

def generate_solinas2010():
    # Load a shor version of the Start.hoc file
    # that does not run the sim
    # it saves the network structure to the Sim_data dir

    if not os.path.exists('SimData'):
        os.mkdir('SimData')
        
    h.load_file('Start_test.hoc')

def compile_mods():
    os.system('nrnivmodl')

    
if __name__ == '__main__':
    compile_mods()
    from neuron import h
    generate_solinas2010()
    structure = loader.load_Solinasetal2010_structure()
