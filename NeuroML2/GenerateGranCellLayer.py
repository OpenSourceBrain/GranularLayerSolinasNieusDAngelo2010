#
#

from neuroml import NeuroMLDocument
from neuroml import Network
from neuroml import Population
from neuroml import Location
from neuroml import Instance
from neuroml import Projection
from neuroml import Connection
from neuroml import IncludeType
from neuroml import InputList
from neuroml import Input
from neuroml import PoissonFiringSynapse

from neuroml import __version__

import neuroml.writers as writers

from pyneuroml import pynml
from pyneuroml.lems.LEMSSimulation import LEMSSimulation

from random import random
from random import seed


def generate_granule_cell_layer(network_id,
                                x_size = 0,     # um
                                y_size = 0,     # um
                                z_size = 0,     # um
                                numCells_grc = 0,
                                numCells_gol = 0,
                                connections = True,
                                connections_method = 'random',
                                connection_probability_grc_gol =   0.2,
                                connection_probability_gol_grc =   0.1,
                                inputs = False,
                                input_firing_rate = 50, # Hz
                                num_inputs_per_grc = 4,
                                validate = True,
                                random_seed = 1234,
                                generate_lems_simulation = False,
                                duration = 500,  # ms
                                dt = 0.05,
                                temperature="32.0 degC"):

    seed(random_seed)

    nml_doc = NeuroMLDocument(id=network_id)

    net = Network(id = network_id, 
                  type = "networkWithTemperature",
                  temperature = temperature)
                  
    net.notes = "Network generated using libNeuroML v%s"%__version__
    nml_doc.networks.append(net)

    # The names of the cell type/component used in the populations (Cell Type in neuroConstruct)
    grc_group_component = "Granule_98"
    gol_group_component = "Golgi_98"

    nml_doc.includes.append(IncludeType(href='%s.cell.nml'%grc_group_component))
    nml_doc.includes.append(IncludeType(href='%s.cell.nml'%gol_group_component))

    # The names of the Exc & Inh groups/populations (Cell Group in neuroConstruct)
    grc_group = "Grans" 
    gol_group = "Golgis" 

    # The names of the network connections 
    net_conn_grc_gol = "NetConn_Grans_Golgis"
    net_conn_gol_grc = "NetConn_Golgis_Grans"

    # The names of the synapse types (should match names at Cell Mechanism/Network tabs in neuroConstruct)
    grc_gol_syn = "AMPA_GranGol"
    gol_grc_syn = "GABAA"

    for syn in [grc_gol_syn, gol_grc_syn]:
        nml_doc.includes.append(IncludeType(href='%s.synapse.nml'%syn))

        
    if network_id == 'Solinas2010':
        # Set size and cell numbers for Solinas 2010 network 
        import sys,os
        NEURON_sim_path = '../NEURON'
        sys.path.append(NEURON_sim_path)
        local_path = os.getcwd()
        os.chdir(NEURON_sim_path)
        print os.getcwd()
        import GenerateSolinas2010 as GS2010
        structure = GS2010.GenerateSolinas2010(generate=True)
        os.chdir(local_path)
        goc_data = structure['Golgis']
        grc_data = structure['Granules']
        grc_pos = grc_data['positions']['data']
        goc_pos = goc_data['positions']['data']
        numCells_grc = grc_pos.shape[0]
        numCells_gol = goc_pos.shape[0]
        
    # Generate excitatory cells 

    grc_pop = Population(id=grc_group, component=grc_group_component, type="populationList", size=numCells_grc)
    net.populations.append(grc_pop)

    for i in range(0, numCells_grc) :
            index = i
            inst = Instance(id=index)
            grc_pop.instances.append(inst)
            inst.location = Location(x=str(grc_pos[i,0]), y=str(grc_pos[i,1]), z=str(grc_pos[i,2]))

    # Generate inhibitory cells
    gol_pop = Population(id=gol_group, component=gol_group_component, type="populationList", size=numCells_gol)
    net.populations.append(gol_pop)

    for i in range(0, numCells_gol) :
            index = i
            inst = Instance(id=index)
            gol_pop.instances.append(inst)
            inst.location = Location(x=str(goc_pos[i,0]), y=str(goc_pos[i,1]), z=str(goc_pos[i,2]))
    
    if connections:

        proj_grc_gol = Projection(id=net_conn_grc_gol, presynaptic_population=grc_group, postsynaptic_population=gol_group, synapse=grc_gol_syn)
        net.projections.append(proj_grc_gol)
        proj_gol_grc = Projection(id=net_conn_gol_grc, presynaptic_population=gol_group, postsynaptic_population=grc_group, synapse=gol_grc_syn)
        net.projections.append(proj_gol_grc)

        count_grc_gol = 0
        count_gol_grc = 0

        # Generate exc -> *  connections

        def add_connection(projection, id, pre_pop, pre_component, pre_cell_id, pre_seg_id, post_pop, post_component, post_cell_id, post_seg_id):

            connection = Connection(id=id, \
                                    pre_cell_id="../%s/%i/%s"%(pre_pop, pre_cell_id, pre_component), \
                                    pre_segment_id=pre_seg_id, \
                                    pre_fraction_along=0.5,
                                    post_cell_id="../%s/%i/%s"%(post_pop, post_cell_id, post_component), \
                                    post_segment_id=post_seg_id,
                                    post_fraction_along=0.5)

            projection.connections.append(connection)

        # Connect Granule cells to Golgi cells
        for i in range(0, numCells_grc):
            # get targets for grc[i] from the structure dict
            # print 'Granule ', i , structure['Granules']['divergence_to_goc']['data'][i][1:]
            for j in structure['Granules']['divergence_to_goc']['data'][i][1:]:
                add_connection(proj_grc_gol, count_grc_gol, grc_group, grc_group_component, i, 0, gol_group, gol_group_component, j, 0)
                count_grc_gol+=1

        # Connect Golgi cells to Granule cells
        for i in range(0, numCells_gol):
            # get targets glomeruli for goc[i] from the lol in structure dict
            # print 'Golgi ', i
            # print structure['Golgis']['divergence_to_glom']['data'][i][1:]
            for k in structure['Golgis']['divergence_to_glom']['data'][i][1:]:
                # get target granule cells for glom[k] from the lol in structure dict
                # print 'Glom ', k
                # print structure['Glomeruli']['divergence_to_grc']['data'][k]
                for j in structure['Glomeruli']['divergence_to_grc']['data'][k][1:]:
                    # print 'Granule ', j
                    add_connection(proj_gol_grc, count_gol_grc, gol_group, gol_group_component, i, 0, grc_group, grc_group_component, j, 0)
                    count_gol_grc+=1

    if inputs:
        
        mf_input_syn = "MF_AMPA"
        nml_doc.includes.append(IncludeType(href='%s.synapse.nml'%mf_input_syn))
        
        rand_spiker_id = "input50Hz"
        
        
        #<poissonFiringSynapse id="Input_8" averageRate="50.0 per_s" synapse="MFSpikeSyn" spikeTarget="./MFSpikeSyn"/>
        pfs = PoissonFiringSynapse(id="input50Hz",
                                   average_rate="%s per_s"%input_firing_rate,
                                   synapse=mf_input_syn,
                                   spike_target="./%s"%mf_input_syn)
                                   
        nml_doc.poisson_firing_synapses.append(pfs)
        
        input_list = InputList(id="Input_0",
                             component=rand_spiker_id,
                             populations=grc_group)
                             
        count = 0
        for i in range(0, numCells_grc):
            
            for j in range(num_inputs_per_grc):
                input = Input(id=count, 
                              target="../%s/%i/%s"%(grc_group, i, grc_group_component), 
                              destination="synapses")  
                input_list.input.append(input)
            
            count += 1
                             
        net.input_lists.append(input_list)


    #######   Write to file  ######    

    print("Saving to file...")
    nml_file = network_id+'.net.nml'
    writers.NeuroMLWriter.write(nml_doc, nml_file)

    print("Written network file to: "+nml_file)


    if validate:

        ###### Validate the NeuroML ######    

        from neuroml.utils import validate_neuroml2
        validate_neuroml2(nml_file) 
        
    if generate_lems_simulation:
        # Create a LEMSSimulation to manage creation of LEMS file
        
        ls = LEMSSimulation("Sim_%s"%network_id, duration, dt)

        # Point to network as target of simulation
        ls.assign_simulation_target(net.id)
        
        # Include generated/existing NeuroML2 files
        ls.include_neuroml2_file('%s.cell.nml'%grc_group_component)
        ls.include_neuroml2_file('%s.cell.nml'%gol_group_component)
        ls.include_neuroml2_file(nml_file)
        

        # Specify Displays and Output Files
        disp_grc = "display_grc"
        ls.create_display(disp_grc, "Voltages Granule cells", "-95", "-38")

        of_grc = 'Volts_file_grc'
        ls.create_output_file(of_grc, "v_grc.dat")
        
        disp_gol = "display_gol"
        ls.create_display(disp_gol, "Voltages Golgi cells", "-95", "-38")

        of_gol = 'Volts_file_gol'
        ls.create_output_file(of_gol, "v_gol.dat")

        for i in range(numCells_grc):
            quantity = "%s/%i/%s/v"%(grc_group, i, grc_group_component)
            ls.add_line_to_display(disp_grc, "GrC %i: Vm"%i, quantity, "1mV", pynml.get_next_hex_color())
            ls.add_column_to_output_file(of_grc, "v_%i"%i, quantity)
            
        for i in range(numCells_gol):
            quantity = "%s/%i/%s/v"%(gol_group, i, gol_group_component)
            ls.add_line_to_display(disp_gol, "Golgi %i: Vm"%i, quantity, "1mV", pynml.get_next_hex_color())
            ls.add_column_to_output_file(of_gol, "v_%i"%i, quantity)

        # Save to LEMS XML file
        lems_file_name = ls.save_to_file()
        
    print "-----------------------------------"


    
if __name__ == "__main__":
    
    
    generate_granule_cell_layer("Solinas2010",
                                connections_method = '../NEURON/',
                                generate_lems_simulation = True)
                                
