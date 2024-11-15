####################################################################################################################################################################
#       Filename:       mainScript.py
#       Summary:        This script simulates a circuit provided as a .net file, outputting the results of the simulation in a .csv file, and possibly graph plots.
#       Description:    Taking in a .csv file it is read and a corresponding virtual cascaded circuit is composed. Its charecteristic two-port ABCD matrix is calc- 
#                       ulated for each frequency defined in the .net file. It is then simulated in tangent to the input and output circuits provided in the .net.
#                       The results of the simulation are stored in a matrix and output in a .csv.
#       Author:         R. Klotins
####################################################################################################################################################################
#                                                               TABLE OF CONTENTS
#
#           Lines     34-45:      Libraries
#
#           Lines    46-177:      Accesibillity functions
#                                 |_ get_n1, get_n2, get_val, get_ABCD.
#
#           Lines   178-276:      Functional Block 1: Cascaded Circuit Compisition                        
#
#           Lines   277-362:      Functional Block 2: Frequency Point Generation
#
#           Lines   363-443:      Functional Block 3: Source and Load Circuits Retrieval
#
#           Lines   444-486:      Functional Block 4: Cascaded Circuit ABCD Matrix Calculations
#
#           Lines   486-596:      Functional Block 5: Composing Bare Output Array
#
#           Lines  596-1651:      Functional Block 6: Calculations and Output Array Population
#                                 |_ calculations, calc_Vin, calc_Iin, calc_Pin, calc_Zin, calc_Vout, calc_Iout, calc_Pout, calc_Zout, calc_Av, calc_Ai
#           
#           Lines 1652-1737:      Main Script
#
#
####################################################################################################################################################################

#=========== Libraries  ============================================================================================================================================

import re
import numpy as np
import math
import cmath
import csv
import sys
import argparse
import matplotlib.pyplot as plt


#=========== Accesiblity functions =================================================================================================================================

def get_n1(component_def):
    '''
    Given a component defintions, returns the integer value of its node 1 (n1).

    param component_def: string component definition (eg. "n1=3 n2=2 R=3")
    
    returns: integer value of n1 (eg. 3). Returns -1 if not found 
    '''
    match = re.search(r'n1\s*=\s*([\d.e+-]+)', component_def)
    if match:
        return int(match.group(1))
    else:
        return -1


def get_n2(component_def):
    '''
    Given a component defintions, returns the integer value of its node 2 (n2).

    param component_def: the string component definition (eg. "n1=1 n2=4 R=3")
    
    returns: integer value of n2 (eg. 4). If not found returns -1.
    '''
   
    match = re.search(r'n2\s*=\s*([\d.e+-]+)', component_def)
    if match:
        return int(match.group(1))
    else:
        return -1
    
def get_val(component_def) :
    '''
    Returns the value and type of the component definition provided. The type can either be
    1: R, 2: G, 3: L, 4: C. If not found, returns -1.

    param component_def: the string component definition (eg. "n1=1 n2=2 C=3")
    
    returns: float value of the compontn in the first position of vector, and the type in the second (eg. [3,4])   

    raises Exception: if values of component cannot be found "No value found for component." 
    '''

    component_type = -1 # If none of the definition types are matched, return -1 (the "not matched flag")
    component_val = -1

    match = re.search(r"R\s*=\s*([\d.e+-]+)", component_def)  # Resistance: component type 1
    if match:
        component_val = float(match.group(1))
        component_type = 1
    match = re.search(r"G\s*=\s*([\d.e+-]+)", component_def)  # Conductance: component type 2
    if match:
        component_val = float(match.group(1))
        component_type = 2
    match = re.search(r"L\s*=\s*([\d.e+-]+)", component_def)  # Inductance: component type 3
    if match:
        component_val = float(match.group(1))
        component_type = 3
    match = re.search(r"C\s*=\s*([\d.e+-]+)", component_def)  # Capacitance: component type 4
    if match:
        component_val = float(match.group(1))
        component_type = 4

    if component_type == -1 or component_val == -1:
        raise Exception("No value found for component")
    
    return [component_val,component_type]   


def get_ABCD(component_def,frequency):
    '''
    Given a component definition and frequency, return the A, B, C, D values of the components ABCD matrix.

    param component_def: String component definition

    param freq: float value of frequency you want to evaluate at

    returns: the 4 matrix values in an array, [A,B,C,D]. 

    raises Expetion: If a component has a zero value, to avoid NaN values and miscalculations "Zero value for components not permitted."     
    '''
    # Use the get_val function to parse the definition and export the component type and value
    component_val = get_val(component_def)
    if component_val[0] == -1 or component_val[1] == -1:    # If either the component type or value couldnt be found, return -1s for "Not found"
        return [-1,-1,-1,-1]
    
    if component_val[0] == 0:
        raise Exception("Zero value for components not permitted. Faulty component:", component_def)
    
    # ABCD Matrix theory:
    #
    # Get the impedance of the value: X
    # If it is a parallel component (connected to common zero node) form the shunt ABCD:
    #[[1,0],[Y,1]]  -> 2 by 2, Y is conductance (1/X)
    # If it is series component (if criteria for parallel not met) form the series ABCD:
    # [[1,X],[0,1]] 
    # 

    # An array of size 4 and complex data-type
    ABCD_out = np.zeros(4,dtype=complex)
    # Resisstance value
    if component_val[1] == 1:
        if (get_n1(component_def) == 0) or (get_n2(component_def) == 0):    # If the component is parallel, one of the nodes will connect to zero, so make shunt ABCD
            ABCD_out = [1,0,(1 / component_val[0]),1]
        else:                                                               # Otherwise make the series impedance 
            ABCD_out = [1,component_val[0],0,1]
    # Conductance value
    if component_val[1] == 2: 
        if (get_n1(component_def) == 0) or (get_n2(component_def) == 0):    
            ABCD_out = [1,0,(component_val[0]),1]                            # resistance = 1 / Conductance
        else:                                                                  
            ABCD_out = [1,(1 / component_val[0]),0,1]                        
    # Inductance value
    if component_val[1] == 3:
        omega = 1j * 2 * cmath.pi * frequency                                # Reactance of inductor eq.
        Z = (omega * component_val[0])
        if (get_n1(component_def) == 0) or (get_n2(component_def) == 0):
            ABCD_out = [1,0,(1 / Z),1]
        else:                                                               
            ABCD_out = [1,(Z),0,1]
    # Capacitance value                                                      # Reactance of capacitor eq.
    if component_val[1] == 4:
        omega = 2 * cmath.pi * frequency
        Z = -1j / (omega * component_val[0])
        if (get_n1(component_def) == 0) or (get_n2(component_def) == 0):
            ABCD_out = [1,0,(1 / Z),1]
        else:                                                               
            ABCD_out = [1,(Z),0,1]
    
    return ABCD_out
        
#=========== Functional Block 1: Cascaded Circuit Composition ======================================================================================================

def cascaded_circuit(file):
    '''
    Takes in the .net file in string array from, every line being in a new element. Retruns the 
    component defintions in a string array which is sorted according to how they would appear 
    in a cascaded circuit diagram (closer to the input -> lower index, vice versa). 
    
    param file: string array, each element is line from the file.

    returns: string array of component definitions, sorted.    
    '''
    
    n_lines = len(file)

    # Looking for the beginning and end of the circuit block, raising errors if
    # they are not found
    for i in range(n_lines):
        if "<CIRCUIT>" in file[i][0:9]: # is "<CIRCUIT>" in line i characters 0-8?
            CIRCUIT_start = i+1         # if it is: the next line is the start of 
            break                       # the terms block
        else:
            CIRCUIT_start = -1
        
    if (CIRCUIT_start == -1):
        raise Exception("Start of CIRCUIT block not found.")

    for i in range(n_lines):
        if "</CIRCUIT>" in file[i][0:10]:
            CIRCUIT_end = i
            break
        else:
            CIRCUIT_end = -1
    
    if (CIRCUIT_end == -1):
        raise Exception("End of CIRCUIT block not found.")

    # Retrieve just the CIRCUIT block
    CIRCUIT_block = file[CIRCUIT_start:CIRCUIT_end]
    
    # Place component definitions in either 'series' or 'parallel' arrays
    series = []
    parallel = []   
    for i in range(len(CIRCUIT_block)):
        # Retrieve the node values from the component in questions, 
        # if that line has no node values or starts with comment then ignore it
        node1 = get_n1(CIRCUIT_block[i])
        node2 = get_n2(CIRCUIT_block[i])
        if (node1 == -1) | (node2 == -1) | (CIRCUIT_block[i][0] == "#"):
            ignore = 1
        else:
            ignore = 0

        if (node1 != 0) & (node2 !=0) & (ignore == 0):          # This part takes advantage of the fact that parallel components will be
            series.append(CIRCUIT_block[i])                     # connected to the zero node to differentiate them from series.
        elif ((node1 == 0) | (node2 == 0)) & (ignore == 0):     
            parallel.append(CIRCUIT_block[i])

    
    # Place the series components in an array, which is the size of the 'parallel' and 'series' arrays combined, starting from the component 
    # connected to node 1, then the component connected to that component, then the one connected to that one etc. This will place the
    # series components in the order that they would appear in on a cascaded diagram. 
    cascaded = [None] * ( len(series) + len(parallel) )
    last_node = 1 
    for i in range( len(series) ):
        for j in range( len(series) ):
            if get_n1(series[j]) == last_node:      # Find the node of the component which is connected to the last node.
                last_node = get_n2(series[j])       # Then set 'last_node' to the other node it is connected to.    
                cascaded[i] = series[j]
                series[j] = "Empty."                # Remove that component definition from the series array.        
                break
            elif get_n2(series[j]) == last_node:    # Same as above but looking at node 2. Because nodes do not have to be
                last_node = get_n1(series[j])       # in order as they appear in the diagram, this allows flexibility for 
                cascaded[i] = series[j]             # the user.
                series[j] = "Empty."
                break

    first_node = 1
    # Inserting the parallel components inbetween the series components.
    for i in range(len(parallel)):
        
        # If the parallel component is connected to the first node, shift 'cascaded' to the right and insert it into the position 0
        if (get_n1(parallel[i]) == first_node) | (get_n2(parallel[i]) == first_node):
            cascaded.insert(0,parallel[i])
            cascaded.pop()

        # Otherwise, use the non-zero node of the parallel component to find where it appears first in the cascaded circuit and place after
        # (A parallel component must be placed inbetween two series components sharing the same node, that means you can place it after the 
        # first component that has the same non-zero node as it)
        else:
            non_zero_node = get_n1(parallel[i]) + get_n2(parallel[i])
            for j in range(len(cascaded)-1):
                if (get_n1(cascaded[j]) == non_zero_node) | (get_n2(cascaded[j]) == non_zero_node): 
                    cascaded.insert(j+1,parallel[i])    # Insert it into position, sliding everything to the right
                    cascaded.pop()                      # Remove (pop) the last index from the array, this is by
                    break                               # necessity an empty cell.

    return cascaded

#=========== Functional Block 2: Frequency Point Generation ======================================================================================================

def frequencies(file):
    '''
    Parse the file and find the TERMS block from which the frequency definitions are extracted.
    From these definitions generate either linear or logarithmically spaced frequency points. 
    
    param file: a string array of .net file, every line is its own element

    returns: a float array of frequency points
    '''
    
    n_lines = len(file)

    # Looking for the beginning and end of the terms block, raising errors if not found
    for i in range(n_lines):
        if "<TERMS>" in file[i][0:8]: # is "<TERMS>" in line i characters 0-8?
            TERMS_start = i+1         # if it is: the next line is the start of 
            break                       # the terms block
        else:
            TERMS_start = -1
        
    if (TERMS_start == -1):
        raise Exception("Start of TERMS block not found.")

    for i in range(n_lines):
        if "</TERMS>" in file[i][0:9]:
            TERMS_end = i
            break
        else:
            TERMS_end = -1
    
    if (TERMS_end == -1):
        raise Exception("End of TERMS block not found.")

    for i in range(n_lines):
        if "</TERMS>" in file[i][0:9]:
            TERMS_end = i
            break
        else:
            TERMS_end = -1

    if (TERMS_end == -1):
        raise Exception("End of TERMS block not found.")
    
    TERMS_block = file[TERMS_start:TERMS_end]

    # Look, line by line, for either the normal or the logarithmic definition of frequency and create the frequency array accordingly.
    for i in range(len(TERMS_block)):

        # Make the regular expressions to find "LFstart", "LFend", "Nfreqs"
        # (we first try catch the log definition to make code simpler. Since
        # the normal definition is contained within the log one, eg. "Fstart"
        # is present in "LFstart".)
        LFstart_regex = r"LFstart\s*=\s*([\d.e+-]+)"
        LFend_regex = r"LFend\s*=\s*([\d.e+-]+)"
        Nfreqs_regex = r"Nfreqs\s*=\s*([\d.e+-]+)"
        # Find out if they match anywhere in the line (file[i])
        LFstart_match = re.search(LFstart_regex, TERMS_block[i])
        LFend_match = re.search(LFend_regex, TERMS_block[i])
        Nfreqs_match = re.search(Nfreqs_regex, TERMS_block[i])
        # If they match then extract values and create the array
        if LFstart_match and LFend_match and Nfreqs_match:
            LFstart = float(LFstart_match.group(1))
            LFend   = float(LFend_match.group(1))
            Nfreqs  =   int(Nfreqs_match.group(1))
                        
            return np.logspace(np.log10(LFstart), np.log10(LFend), num=Nfreqs)


        #  If they do not match try "Fstart", "Fend", and "Nfreq"
        Fstart_regex = r"Fstart\s*=\s*([\d.e+-]+)"
        Fend_regex = r"Fend\s*=\s*([\d.e+-]+)"
        Nfreqs_regex = r"Nfreqs\s*=\s*([\d.e+-]+)"
        # Find out if they match anywhere in the line (file[i])
        Fstart_match = re.search(Fstart_regex, TERMS_block[i])
        Fend_match = re.search(Fend_regex, TERMS_block[i])
        Nfreqs_match = re.search(Nfreqs_regex, TERMS_block[i])
        # If they match then extract values and create the array
        if Fstart_match and Fend_match and Nfreqs_match:
            Fstart  = float(Fstart_match.group(1))
            Fend    = float(Fend_match.group(1)  )
            Nfreqs  =   int(Nfreqs_match.group(1))
            return np.linspace(Fstart, Fend, Nfreqs)
        

#=========== Functional Block 3: Source and Load Circuits Retrieval ======================================================================================================7


def get_VsRsRl(file):
    '''
    From the <TERMS> block find the definitions of the load and input circuits, operate on them
    if necessary to normalise (turning GS into RS). Then return VS, Rs, Rl.

    param file: a string array of .net file, each element is a new line

    returns: float array of VT, RS, RL values

    rasies Exception: if the terms cannot be found "Terms values not present or not formatted correctly."
                      or if the <TERMS> block is not found "Start/End of TERMS block not found.") 
    '''
    n_lines = len(file)

    # Looking for the beginning and end of the terms block, raising errors if not found
    for i in range(n_lines):
        if "<TERMS>" in file[i][0:8]: # is "<TERMS>" in line i characters 0-8?
            TERMS_start = i+1         # if it is: the next line is the start of 
            break                       # the terms block
        else:
            TERMS_start = -1 
    if (TERMS_start == -1):
        raise Exception("Start of TERMS block not found.")

    for i in range(n_lines):
        if "</TERMS>" in file[i][0:9]:
            TERMS_end = i
            break
        else:
            TERMS_end = -1
    if (TERMS_end == -1):
        raise Exception("End of TERMS block not found.")

    TERMS_block = file[TERMS_start:TERMS_end]
    # The values being -1 means they are not found yet.
    RS = -1
    RL = -1
    IN = -1
    VT = -1

    for i in range(len(TERMS_block)):
        line = TERMS_block[i]
        # Ignore line if it starts with "#"
        if line[0] == "#":
            continue

        # If RS found: RS = val
        match = re.search(r'RS\s*=\s*([\d.e+-]+)', line)
        if match:
            RS = float(match.group(1))
        # If GS found: RS = 1/val
        match = re.search(r'GS\s*=\s*([\d.e+-]+)', line)
        if match:
            RS = 1/float(match.group(1))
        # If IN found: IN = val, VT = -1
        match = re.search(r'IN\s*=\s*([\d.e+-]+)', line)
        if match:
            IN = float(match.group(1))
        # If VT found: VT = val, IN = -1
        match = re.search(r'VT\s*=\s*([\d.e+-]+)', line)
        if match:
            VT = float(match.group(1))
        # If RL found: RL = val
        match = re.search(r'RL\s*=\s*([\d.e+-]+)', line)
        if match:
            RL = float(match.group(1))

    # If the neccessary values have not been found then raise exception. (either: RS or RL or (both VT and IN )
    if (RS == -1) or (RL == -1) or ((IN == -1) and (VT == -1)):
        raise Exception("Terms values not present or not formatted correctly.")
    
    # Check if VT was not found, meaning that a Norton source was used, then convert IN to Thevanin
    if VT == -1:
        VT = IN * RS

    return [VT,RS,RL]

       
#=========== Functional Block 4: Cascaded Circuit ABCD Matrix Calculations ======================================================================================================7
   
def transfer_functions(cascaded, frequencies):
    '''
    For every frequency, go through the cascaded circuit and multiply together each ABCD matrix till the end to get 
    the total transfer fucntion of the circuit, and reutrn it.
    
    param cascaded: String array of componenet definitions ordered according to how a cascaded circuit would be.

    param frequencies: Float array of the frequency points desired by the user.

    returns: The float A, B, C, D values of each frequency in 4 one dimentional arrays. 
    '''
    n_frequencies = len(frequencies)

    # Make a complex array for A,B,C,D values of the ABCD matrix.
    A_f = np.zeros(n_frequencies, dtype=complex)
    B_f = np.zeros(n_frequencies, dtype=complex)
    C_f = np.zeros(n_frequencies, dtype=complex)
    D_f = np.zeros(n_frequencies, dtype=complex)
    # For every frequency, calculate total ABCD (T)
    for i in range(n_frequencies):
        for j in range(len(cascaded)):
            if j == 0:
                # For the first component, initialise the transfer function, which you will multiply by the rest to get the total ABCD matrix.
                [A,B,C,D] = get_ABCD(cascaded[0],frequencies[i])
                T = np.array([[A, B], [C, D]])
                continue
            #If we are not looking at the first component:
            # Find the ABCD matrix of the current component and multiply by the running current total variable (T) 
            [A,B,C,D] = get_ABCD(cascaded[j],frequencies[i])
            ABCD = np.array([[A,B],[C,D]])
            T = np.dot(T, ABCD) # Dot multiplication of the total transfer function so far (T) and the current components ABCD matrix
        # Save this particular frequencies (frequency i) transfer function ABCD matric values.
        A_f[i] = T[0][0]
        B_f[i] = T[0][1]
        C_f[i] = T[1][0]
        D_f[i] = T[1][1] 

    return A_f,B_f,C_f,D_f


#=========== Functional Block 5: Composing Bare Output Array =============================================================================================================================================

def output_array(file, frequencies):
    '''
    Extraxt the desired output parameters from the file, place these in an array for later use. Make a 2d array with a length of one greater than the number 
    of values in the desired output file, and a height of two greater than the number of frequencies. Populate first column with frequencies and their labels.

    param file: The input file as a string array, each element being a seperate line

    param frequencies: float array holding the frequency values

    return: An array of arrays, each sub-array is a row, dimensions according to number of parameters and frequencies. And a string array containging the def-
            initions of the desired parameters (eg "Vin V") in positions 2n+1 and their labels (eg. "Vin") in positions 2n.

    raises Exceptions: If the start or end of the <OUTPUT> block not found "Start/End of OUTPUT block not found."
    '''

    n_frequencies = len(frequencies)
    n_lines = len(file)

    # Looking for the beginning and end of the terms block, raising errors if not found
    for i in range(n_lines):
        if "<OUTPUT>" in file[i][0:8]: # is "<TERMS>" in line i characters 0-8?
            OUTPUT_start = i+1         # if it is: the next line is the start of 
            break                       # the terms block
        else:
            OUTPUT_start = -1 
    if (OUTPUT_start == -1):
        raise Exception("Start of OUTPUT block not found.")

    for i in range(n_lines):
        if "</OUTPUT>" in file[i][0:9]:
            OUTPUT_end = i
            break
        else:
            OUTPUT_end = -1
    if (OUTPUT_end == -1):
        raise Exception("End of OUTPUT block not found.")

    OUTPUT_block = file[OUTPUT_start:OUTPUT_end]

    # Look through the <OUTPUT> block and finds all the output parameters requested.
    # It then places the parameters label (eg. 'Vin') into the parameters array then
    # it places the line it was found from into the parameters array. This is so that
    # if there are multiple definitions on one line, they can be still be found, one
    # will not shadow the other. So the number of parameters will be the length of the 
    # array divided by two.
    parameters = []
    for i in range(len(OUTPUT_block)):
        # If the line starts with a hash then ignore.
        if  OUTPUT_block[i][0] == "#":
            continue
        match = re.search(r'Vin',  OUTPUT_block[i])     
        if match:
            parameters.append('Vin')
            parameters.append(OUTPUT_block[i])
        match = re.search(r'Iin',  OUTPUT_block[i])     
        if match:
            parameters.append('Iin')
            parameters.append(OUTPUT_block[i])
        match = re.search(r'Pin',  OUTPUT_block[i])     
        if match:
            parameters.append('Pin')
            parameters.append(OUTPUT_block[i])
        match = re.search(r'Zin',  OUTPUT_block[i])     
        if match:
            parameters.append('Zin')
            parameters.append(OUTPUT_block[i])
        match = re.search(r'Vout',  OUTPUT_block[i])     
        if match:
            parameters.append('Vout')
            parameters.append(OUTPUT_block[i])
        match = re.search(r'Iout',  OUTPUT_block[i])     
        if match:
            parameters.append('Iout')
            parameters.append(OUTPUT_block[i])
        match = re.search(r'Pout',  OUTPUT_block[i])     
        if match:
            parameters.append('Pout')
            parameters.append(OUTPUT_block[i])
        match = re.search(r'Zout',  OUTPUT_block[i])     
        if match:
            parameters.append('Zout')
            parameters.append(OUTPUT_block[i])
        match = re.search(r'Av',  OUTPUT_block[i])     
        if match:
            parameters.append('Av')
            parameters.append(OUTPUT_block[i])
        match = re.search(r'Ai',  OUTPUT_block[i])     
        if match:
            parameters.append('Ai')
            parameters.append(OUTPUT_block[i])

    n_parameters = len(parameters)

    sig_figs = 4    # Formats the frequencies to be to 4 significant figures and in exponent form.
    frequencies = [f"{x:.{sig_figs-1}e}" for x in frequencies]

    # Rows: number of frequencies + 2
    # Columns: number of parameters + 1
    array_out = [[0 for i in range(n_parameters+1)] for j in range(n_frequencies+2)]

    array_out[0][0] =  ("Freq").rjust(10)
    array_out[1][0] = ("Hz").rjust(10)
    for i in range(n_frequencies):
        array_out[i+2][0] = str(frequencies[i]).rjust(10)

    return array_out, parameters


#=========== Functional Block 6: Calculations and Output Array Population ===============================================================================================================================================

def calculations(A,B,C,D, output_array, parameter_list, VsRsRl):
    '''
    Takes in a list of parameters that need to be calculated, with their definition functions, aswell as the circuits
    ABCD matricies and frequencies. Uses these to calculate the values of the parameters while filling in the output
    array. The filled out and labelled output array is returned. The calculation functions are below this one.

    param A,B,C,D: 4 float arrays containing the imaginary values of the cascaded circuits ABCD matricies for each frequency

    param output_array: An array of arrays, with first column filled out

    param parameter_list: A string array of parameter labels and definitions telling the program which values to calculate (eg ["Vin", "Vin dBV"])

    returns: An array of arrays, fully filled out with required values. Each one dimensional array is a row in table. 

    raises Exception: If 2nth position doesn't hold a parameter label "The array containing parameters is not correctly formulated."
    '''

    n_parameters = int(len(parameter_list) / 2)
    
    for i in range(n_parameters):
        # The [2 x nth] value is the label of the parameter, the [2 x nth + 1] is the parameters definition.
        # Here we see what the label is and call the corresponding 'get_Xxx' function, which returns the two columns, their labels and units.
        # { Clarity: 2*i + 1 will be the parameter definition whereas 2*i is the parameter flag (flag: "Vin", def: "Vin dBmV") }
        if parameter_list[2*i] == "Vin":
            left_column,right_column,labels,units = calc_Vin(VsRsRl, A,B,C,D, parameter_list[2*i+1]) 
        elif parameter_list[2*i] == "Iin":              
            left_column,right_column,labels,units = calc_Iin(VsRsRl, A,B,C,D, parameter_list[2*i+1])
        elif parameter_list[2*i] == "Pin":
            left_column,right_column,labels,units = calc_Pin(VsRsRl, A,B,C,D, parameter_list[2*i+1])
        elif parameter_list[2*i] == "Zin":
            left_column,right_column,labels,units = calc_Zin(VsRsRl, A,B,C,D, parameter_list[2*i+1])
        elif parameter_list[2*i] == "Vout":
            left_column,right_column,labels,units = calc_Vout(VsRsRl, A,B,C,D, parameter_list[2*i+1])
        elif parameter_list[2*i] == "Iout":
            left_column,right_column,labels,units = calc_Iout(VsRsRl, A,B,C,D, parameter_list[2*i+1])
        elif parameter_list[2*i] == "Pout":
            left_column,right_column,labels,units = calc_Pout(VsRsRl, A,B,C,D, parameter_list[2*i+1])
        elif parameter_list[2*i] == "Zout":
            left_column,right_column,labels,units = calc_Zout(VsRsRl, A,B,C,D, parameter_list[2*i+1])
        elif parameter_list[2*i] == "Av":
            left_column,right_column,labels,units = calc_Av(VsRsRl, A,B,C,D, parameter_list[2*i+1])
        elif parameter_list[2*i] == "Ai":
            left_column,right_column,labels,units = calc_Ai(VsRsRl, A,B,C,D, parameter_list[2*i+1])
        else: # If none of the parameters are present in the '2*nth' position, it is incorrectly formulated. 
            raise Exception("The array containing parameters is not correctly formulated.") 
    
        # Now that you have the two columns containing the parameters values and labels, place them
        # into the next available spaces in the output array.
        for i in range(len(output_array[0])): # For the width of the output array
            if output_array[0][i] == 0:
                empty_space = i
                break
        
        # Place labels and units at the top of the 2 new empty columns
        output_array[0][empty_space]    = labels[0].rjust(11)
        output_array[0][empty_space+1]  = labels[1].rjust(11)
        output_array[1][empty_space]    = units[0].rjust(11)
        output_array[1][empty_space+1]  = units[1].rjust(11)
        
        # Populate left column
        for i in range(len(left_column)):
            output_array[2+i][empty_space] = str(left_column[i]).rjust(11)
        # Populate right column
        for i in range(len(right_column)):
            output_array[2+i][empty_space+1] = str(right_column[i]).rjust(11)    # Here I right justify the values and labels so the columns align

    return output_array

#======================== Calculation functions for individual parameters ========================================================================================================================================
def calc_Vin(VsRsRl, A,B,C,D, parameter_def):
    '''
    Using the transfer function matricies of the cascaded circuit as well as source and load circuits calculate the values
    for Vin. These are formatted and output as 2 columns.

    param VsRsRl: Array containing [Source voltage, source resistance, load resistance]

    param A,B,C,D: 4 arrays containing the corresponding values for the transfer matricies of each frequency.

    param parameter_def: the definition of the Vin parameter, indicating how it ought to be formatted

    returns: two columns (2 by n_frequencies array)
    '''
    Vs = VsRsRl[0]  # Thevanin source, Thevanin resistance, Load resistance
    Rs = VsRsRl[1]
    Rl = VsRsRl[2]
    n_freq = len(A)
    # Create two empty columns for real and imaginary values, will format and scale later.
    Re_Vin = np.zeros(n_freq,dtype = complex)
    Im_Vin = np.zeros(n_freq,dtype = complex)
    # Create the labels for the parameter and unit
    label_parameter = "Vin"
    label_unit      = "V"
    # For every fr@equency calculate the Vin value
    for i in range(n_freq):
        Zin = (A[i]*Rl + B[i]) / (C[i]*Rl + D[i])
        Vin = (Zin*Vs)/(Zin + Rs)
        Re_Vin[i] = Vin.real
        Im_Vin[i] = Vin.imag

    # Now check for prefixes and scale
    match = re.search(r'pV', parameter_def)
    if match:
        Re_Vin = Re_Vin * 10**12
        Im_Vin = Im_Vin * 10**12
        label_unit = "p" + label_unit
    match = re.search(r'nV', parameter_def)
    if match:
        Re_Vin = Re_Vin * 10**9
        Im_Vin = Im_Vin * 10**9
        label_unit = "n" + label_unit
    match = re.search(r'uV', parameter_def)
    if match:
        Re_Vin = Re_Vin * 10**6
        Im_Vin = Im_Vin * 10**6
        label_unit = "u" + label_unit
    match = re.search(r'mV', parameter_def)
    if match:
        Re_Vin = Re_Vin * 10**3
        Im_Vin = Im_Vin * 10**3
        label_unit = "m" + label_unit
    match = re.search(r'kV', parameter_def)
    if match:
        Re_Vin = Re_Vin * 10**-3
        Im_Vin = Im_Vin * 10**-3
        label_unit = "k" + label_unit
    match = re.search(r'MV', parameter_def)
    if match:
        Re_Vin = Re_Vin * 10**-6
        Im_Vin = Im_Vin * 10**-6
        label_unit = "M" + label_unit
    match = re.search(r'GV', parameter_def)
    if match:
        Re_Vin = Re_Vin * 10**-9
        Im_Vin = Im_Vin * 10**-9
        label_unit = "G" + label_unit

    # If there is a dB prefix before the V unit, convert to the dB format and amend labels,
    match = re.search(r'dB.?V', parameter_def)
    if match:
        complex_values = Re_Vin + Im_Vin * 1j   
        out_left = 20 * np.log10( np.abs(complex_values) )       # Populate output columns correctly. (to find dB: 10log10( abs of complex value ) )
        [max(x, -160) for x in out_left]                            # Limiting dB output to be greater than -160
        out_right = np.angle(complex_values)
        label_parameter_left    = "|" + label_parameter + "|"       # Compose correct labels for the dB case
        label_parameter_right   = "/_" + label_parameter
        label_unit_left         = "dB" + label_unit
        label_unit_right        = "Rads"
    else:
        out_left = Re_Vin       # Output columns
        out_right = Im_Vin
        label_parameter_left    = "Re(" + label_parameter + ")"        # Compose correct labels for normal case
        label_parameter_right   = "Im(" + label_parameter + ")"
        label_unit_left         = label_unit
        label_unit_right        = label_unit
    
    labels = [label_parameter_left,label_parameter_right]   # Combine the labels and units into two arrays for simplicity
    units  = [label_unit_left,label_unit_right]
    out_left = out_left.real                                # Return just the real parts of the values
    out_right = out_right.real
    
    sig_figs = 4    # Formats the values to be to 4 significant figures and in exponent form.
    out_left= [f"{x:.{sig_figs-1}e}" for x in out_left]
    out_right = [f"{x:.{sig_figs-1}e}" for x in out_right]
    return out_left,out_right,labels,units

def calc_Iin(VsRsRl, A,B,C,D, parameter_def):
    '''
    Using the transfer function matricies of the cascaded circuit as well as source and load circuits calculate the values
    for Iin. These are formatted and output as 2 columns.

    param VsRsRl: Array containing [Source voltage, source resistance, load resistance]

    param A,B,C,D: 4 arrays containing the corresponding values for the transfer matricies of each frequency.

    param parameter_def: the definition of the Iin parameter, indicating how it ought to be formatted

    returns: two columns (2 by n_frequencies array), and returns the parameter labels [left,right]
             and the unit labels [left,right]
    '''

    Vs = VsRsRl[0]# Thevanin source, Thevanin resistance, Load resistance
    Rs = VsRsRl[1]
    Rl = VsRsRl[2]
    n_freq = len(A)
    # Create two empty columns for real and imaginary values, will format and scale later.
    Re_Iin = np.zeros(n_freq,dtype = complex)
    Im_Iin = np.zeros(n_freq,dtype = complex)
    # Create the labels for the parameter and unit
    label_parameter = "Iin"
    label_unit      = "A"
    # For every fr@equency calculate the Iin value
    for i in range(n_freq):
        Zin = (A[i]*Rl + B[i]) / (C[i]*Rl + D[i])
        Iin = (Vs)/(Zin + Rs)
        Re_Iin[i] = Iin.real
        Im_Iin[i] = Iin.imag

    # Now check for prefixes and scale
    match = re.search(r'pA', parameter_def)
    if match:
        Re_Iin = Re_Iin * 10**12
        Im_Iin = Im_Iin * 10**12
        label_unit = "p" + label_unit
    match = re.search(r'nA', parameter_def)
    if match:
        Re_Iin = Re_Iin * 10**9
        Im_Iin = Im_Iin * 10**9
        label_unit = "n" + label_unit
    match = re.search(r'uA', parameter_def)
    if match:
        Re_Iin = Re_Iin * 10**6
        Im_Iin = Im_Iin * 10**6
        label_unit = "u" + label_unit
    match = re.search(r'mA', parameter_def)
    if match:
        Re_Iin = Re_Iin * 10**3
        Im_Iin = Im_Iin * 10**3
        label_unit = "m" + label_unit
    match = re.search(r'kA', parameter_def)
    if match:
        Re_Iin = Re_Iin * 10**-3
        Im_Iin = Im_Iin * 10**-3
        label_unit = "k" + label_unit
    match = re.search(r'MA', parameter_def)
    if match:
        Re_Iin = Re_Iin * 10**-6
        Im_Iin = Im_Iin * 10**-6
        label_unit = "M" + label_unit
    match = re.search(r'GA', parameter_def)
    if match:
        Re_Iin = Re_Iin * 10**-9
        Im_Iin = Im_Iin * 10**-9
        label_unit = "G" + label_unit

    # If there is a dB prefix before the A unit, convert to the dB format and amend labels,
    match = re.search(r'dB.?A', parameter_def)
    if match:
        complex_values = Re_Iin + Im_Iin * 1j   
        out_left = 20 * np.log10( np.abs(complex_values) )       # Populate output columns correctly. (to find dB: 10log10( abs of complex value ) )
        [max(x, -160) for x in out_left]                            # Limiting dB output to be greater than -160
        out_right = np.angle(complex_values)
        label_parameter_left    = "|" + label_parameter + "|"       # Compose correct labels for the dB case
        label_parameter_right   = "/_" + label_parameter
        label_unit_left         = "dB" + label_unit
        label_unit_right        = "Rads"
    else:
        out_left = Re_Iin       # Output columns
        out_right = Im_Iin
        label_parameter_left    = "Re(" + label_parameter + ")"        # Compose correct labels for normal case
        label_parameter_right   = "Im(" + label_parameter + ")"
        label_unit_left         = label_unit
        label_unit_right        = label_unit
    
    labels = [label_parameter_left,label_parameter_right]   # Combine the labels and units into two arrays for simplicity
    units  = [label_unit_left,label_unit_right]
    out_left = out_left.real                                # Return just the real parts of the values
    out_right = out_right.real
     
    sig_figs = 4    # Formats the values to be to 4 significant figures and in exponent form.
    out_left= [f"{x:.{sig_figs-1}e}" for x in out_left]
    out_right = [f"{x:.{sig_figs-1}e}" for x in out_right]

    return out_left,out_right,labels,units

def calc_Pin(VsRsRl, A,B,C,D, parameter_def):
    '''
    Using the transfer function matricies of the cascaded circuit as well as source and load circuits calculate the values
    for Pin. These are formatted and output as 2 columns.

    param VsRsRl: Array containing [Source voltage, source resistance, load resistance]

    param A,B,C,D: 4 arrays containing the corresponding values for the transfer matricies of each frequency.

    param parameter_def: the definition of the Pin parameter, indicating how it ought to be formatted

    returns: two columns (2 by n_frequencies array), and returns the parameter labels [left,right]
             and the unit labels [left,right]
    '''

    Vs = VsRsRl[0]# Thevanin source, Thevanin resistance, Load resistance
    Rs = VsRsRl[1]
    Rl = VsRsRl[2]
    n_freq = len(A)
    # Create two empty columns for real and imaginary values, will format and scale later.
    Re_Pin = np.zeros(n_freq,dtype = complex)
    Im_Pin = np.zeros(n_freq,dtype = complex)
    # Create the labels for the parameter and unit
    label_parameter = "Pin"
    label_unit      = "W"
    # For every fr@equency calculate the Pin value (Voltage x Conjugate of Current)
    for i in range(n_freq):        
        Zin = (A[i]*Rl + B[i]) / (C[i]*Rl + D[i])
        Iin = (Vs)/(Zin + Rs)
        Vin = (Zin*Vs)/(Zin + Rs)
        Pin = Vin * np.conj(Iin)
        Re_Pin[i] = Pin.real
        Im_Pin[i] = Pin.imag

    # Now check for prefixes and scale
    match = re.search(r'pW', parameter_def)
    if match:
        Re_Pin = Re_Pin * 10**12
        Im_Pin = Im_Pin * 10**12
        label_unit = "p" + label_unit
    match = re.search(r'nW', parameter_def)
    if match:
        Re_Pin = Re_Pin * 10**9
        Im_Pin = Im_Pin * 10**9
        label_unit = "n" + label_unit
    match = re.search(r'uW', parameter_def)
    if match:
        Re_Pin = Re_Pin * 10**6
        Im_Pin = Im_Pin * 10**6
        label_unit = "u" + label_unit
    match = re.search(r'mW', parameter_def)
    if match:
        Re_Pin = Re_Pin * 10**3
        Im_Pin = Im_Pin * 10**3
        label_unit = "m" + label_unit
    match = re.search(r'kW', parameter_def)
    if match:
        Re_Pin = Re_Pin * 10**-3
        Im_Pin = Im_Pin * 10**-3
        label_unit = "k" + label_unit
    match = re.search(r'MW', parameter_def)
    if match:
        Re_Pin = Re_Pin * 10**-6
        Im_Pin = Im_Pin * 10**-6
        label_unit = "M" + label_unit
    match = re.search(r'GW', parameter_def)
    if match:
        Re_Pin = Re_Pin * 10**-9
        Im_Pin = Im_Pin * 10**-9
        label_unit = "G" + label_unit

    # If there is a dB prefix before the W unit, convert to the dB format and amend labels,
    match = re.search(r'dB.?W', parameter_def)
    if match:
        complex_values = Re_Pin + Im_Pin * 1j   
        out_left = 10 * np.log10( np.abs(complex_values) )       # Populate output columns correctly. (to find dB: 10log10( abs of complex value ) )
        [max(x, -100) for x in out_left]                            # Limiting dB output to be greater than -100
        out_right = np.angle(complex_values)
        label_parameter_left    = "|" + label_parameter + "|"       # Compose correct labels for the dB case
        label_parameter_right   = "/_" + label_parameter
        label_unit_left         = "dB" + label_unit
        label_unit_right        = "Rads"
    else:
        out_left = Re_Pin       # Output columns
        out_right = Im_Pin
        label_parameter_left    = "Re(" + label_parameter + ")"        # Compose correct labels for normal case
        label_parameter_right   = "Im(" + label_parameter + ")"
        label_unit_left         = label_unit
        label_unit_right        = label_unit
    
    labels = [label_parameter_left,label_parameter_right]   # Combine the labels and units into two arrays for simplicity
    units  = [label_unit_left,label_unit_right]
    out_left = out_left.real                                # Return just the real parts of the values
    out_right = out_right.real
    sig_figs = 4    # Formats the values to be to 4 significant figures and in exponent form.
    out_left= [f"{x:.{sig_figs-1}e}" for x in out_left]
    out_right = [f"{x:.{sig_figs-1}e}" for x in out_right]
    
    return out_left,out_right,labels,units

def calc_Zin(VsRsRl, A,B,C,D, parameter_def):
    '''
    Using the transfer function matricies of the cascaded circuit as well as source and load circuits to calculate the values
    for Zin. These are formatted and output as 2 columns.

    param VsRsRl: Array containing [Source voltage, source resistance, load resistance]

    param A,B,C,D: 4 arrays containing the corresponding values for the transfer matricies of each frequency.

    param parameter_def: the definition of the Zin parameter, indicating how it ought to be formatted

    returns: two columns (2 by n_frequencies array), and returns the parameter labels [left,right]
             and the unit labels [left,right]
    '''

    Vs = VsRsRl[0]# Thevanin source, Thevanin resistance, Load resistance
    Rs = VsRsRl[1]
    Rl = VsRsRl[2]
    n_freq = len(A)
    # Create two empty columns for real and imaginary values, will format and scale later.
    Re_Zin = np.zeros(n_freq,dtype = complex)
    Im_Zin = np.zeros(n_freq,dtype = complex)
    # Create the labels for the parameter and unit
    label_parameter = "Zin"
    label_unit      = "Ohms"
    # For every freequency calculate the Zin value (A*Zl + B) / (C*Zl + D)
    for i in range(n_freq):        
        Zin = (A[i]*Rl + B[i]) / (C[i]*Rl + D[i])
        Re_Zin[i] = Zin.real
        Im_Zin[i] = Zin.imag

    # Now check for prefixes and scale
    match = re.search(r'pOhms', parameter_def)
    if match:
        Re_Zin = Re_Zin * 10**12
        Im_Zin = Im_Zin * 10**12
        label_unit = "p" + label_unit
    match = re.search(r'nOhms', parameter_def)
    if match:
        Re_Zin = Re_Zin * 10**9
        Im_Zin = Im_Zin * 10**9
        label_unit = "n" + label_unit
    match = re.search(r'uOhms', parameter_def)
    if match:
        Re_Zin = Re_Zin * 10**6
        Im_Zin = Im_Zin * 10**6
        label_unit = "u" + label_unit
    match = re.search(r'mOhms', parameter_def)
    if match:
        Re_Zin = Re_Zin * 10**3
        Im_Zin = Im_Zin * 10**3
        label_unit = "m" + label_unit
    match = re.search(r'kOhms', parameter_def)
    if match:
        Re_Zin = Re_Zin * 10**-3
        Im_Zin = Im_Zin * 10**-3
        label_unit = "k" + label_unit
    match = re.search(r'MOhms', parameter_def)
    if match:
        Re_Zin = Re_Zin * 10**-6
        Im_Zin = Im_Zin * 10**-6
        label_unit = "M" + label_unit
    match = re.search(r'GOhms', parameter_def)
    if match:
        Re_Zin = Re_Zin * 10**-9
        Im_Zin = Im_Zin * 10**-9
        label_unit = "G" + label_unit

    # If there is a dB prefix before the Ohms unit, convert to the dB format and amend labels,
    match = re.search(r'dB.?Ohms', parameter_def)
    if match:
        complex_values = Re_Zin + Im_Zin * 1j   
        out_left = 20 * np.log10( np.abs(complex_values) )       # Populate output columns correctly. (to find dB: 10log10( abs of complex value ) )
        [max(x, -160) for x in out_left]                            # Limiting dB output to be greater than -160
        out_right = np.angle(complex_values)
        label_parameter_left    = "|" + label_parameter + "|"       # Compose correct labels for the dB case
        label_parameter_right   = "/_" + label_parameter
        label_unit_left         = "dB" + label_unit
        label_unit_right        = "Rads"
    else:
        out_left = Re_Zin       # Output columns
        out_right = Im_Zin
        label_parameter_left    = "Re(" + label_parameter + ")"        # Compose correct labels for normal case
        label_parameter_right   = "Im(" + label_parameter + ")"
        label_unit_left         = label_unit
        label_unit_right        = label_unit
    
    labels = [label_parameter_left,label_parameter_right]   # Combine the labels and units into two arrays for simplicity
    units  = [label_unit_left,label_unit_right]
    out_left = out_left.real                                # Return just the real parts of the values
    out_right = out_right.real
    
    sig_figs = 4    # Formats the values to be to 4 significant figures and in exponent form.
    out_left= [f"{x:.{sig_figs-1}e}" for x in out_left]
    out_right = [f"{x:.{sig_figs-1}e}" for x in out_right]

    return out_left,out_right,labels,units


def calc_Vout(VsRsRl, A,B,C,D, parameter_def):
    '''
    Using the transfer function matricies of the cascaded circuit as well as source and load circuits to calculate the values
    for Vout. These are formatted and output as 2 columns.

    param VsRsRl: Array containing [Source voltage, source resistance, load resistance]

    param A,B,C,D: 4 arrays containing the corresponding values for the transfer matricies of each frequency.

    param parameter_def: the definition of the Vout parameter, indicating how it ought to be formatted

    returns: two columns (2 by n_frequencies array), and returns the parameter labels [left,right]
             and the unit labels [left,right]
    '''

    Vs = VsRsRl[0]# Thevanin source, Thevanin resistance, Load resistance
    Rs = VsRsRl[1]
    Rl = VsRsRl[2]
    n_freq = len(A)
    # Create two empty columns for real and imaginary values, will format and scale later.
    Re_Vout = np.zeros(n_freq,dtype = complex)
    Im_Vout = np.zeros(n_freq,dtype = complex)
    # Create the labels for the parameter and unit
    label_parameter = "Vout"
    label_unit      = "V"
    # For every freequency calculate the Vout value (Vin / A + V Yl) (Yl = 1/Rl)
    for i in range(n_freq):        
        Zin = (A[i]*Rl + B[i]) / (C[i]*Rl + D[i])
        Vin = (Zin*Vs)/(Zin + Rs)
        Vout = Vin / (A[i] + (B[i] / Rl))
        Re_Vout[i] = Vout.real
        Im_Vout[i] = Vout.imag

    # Now check for prefixes and scale
    match = re.search(r'pV', parameter_def)
    if match:
        Re_Vout = Re_Vout * 10**12
        Im_Vout = Im_Vout * 10**12
        label_unit = "p" + label_unit
    match = re.search(r'nV', parameter_def)
    if match:
        Re_Vout = Re_Vout * 10**9
        Im_Vout = Im_Vout * 10**9
        label_unit = "n" + label_unit
    match = re.search(r'uV', parameter_def)
    if match:
        Re_Vout = Re_Vout * 10**6
        Im_Vout = Im_Vout * 10**6
        label_unit = "u" + label_unit
    match = re.search(r'mV', parameter_def)
    if match:
        Re_Vout = Re_Vout * 10**3
        Im_Vout = Im_Vout * 10**3
        label_unit = "m" + label_unit
    match = re.search(r'kV', parameter_def)
    if match:
        Re_Vout = Re_Vout * 10**-3
        Im_Vout = Im_Vout * 10**-3
        label_unit = "k" + label_unit
    match = re.search(r'MV', parameter_def)
    if match:
        Re_Vout = Re_Vout * 10**-6
        Im_Vout = Im_Vout * 10**-6
        label_unit = "M" + label_unit
    match = re.search(r'GV', parameter_def)
    if match:
        Re_Vout = Re_Vout * 10**-9
        Im_Vout = Im_Vout * 10**-9
        label_unit = "G" + label_unit

    # If there is a dB prefix before the V unit, convert to the dB format and amend labels,
    match = re.search(r'dB.?V', parameter_def)
    if match:
        complex_values = Re_Vout + Im_Vout * 1j   
        out_left = 20 * np.log10( np.abs(complex_values) )       # Populate output columns correctly. (to find dB: 10log10( abs of complex value ) )
        [max(x, -160) for x in out_left]                            # Limiting dB output to be greater than -160
        out_right = np.angle(complex_values)
        label_parameter_left    = "|" + label_parameter + "|"       # Compose correct labels for the dB case
        label_parameter_right   = "/_" + label_parameter
        label_unit_left         = "dB" + label_unit
        label_unit_right        = "Rads"
    else:
        out_left = Re_Vout       # Output columns
        out_right = Im_Vout
        label_parameter_left    = "Re(" + label_parameter + ")"        # Compose correct labels for normal case
        label_parameter_right   = "Im(" + label_parameter + ")"
        label_unit_left         = label_unit
        label_unit_right        = label_unit
    
    labels = [label_parameter_left,label_parameter_right]   # Combine the labels and units into two arrays for simplicity
    units  = [label_unit_left,label_unit_right]
    out_left = out_left.real                                # Return just the real parts of the values
    out_right = out_right.real
    sig_figs = 4    # Formats the values to be to 4 significant figures and in exponent form.
    out_left= [f"{x:.{sig_figs-1}e}" for x in out_left]
    out_right = [f"{x:.{sig_figs-1}e}" for x in out_right]
    

    return out_left,out_right,labels,units

def calc_Iout(VsRsRl, A,B,C,D, parameter_def):
    '''
    Using the transfer function matricies of the cascaded circuit as well as source and load circuits to calculate the values
    for Iout. These are formatted and output as 2 columns.

    param VsRsRl: Array containing [Source voltage, source resistance, load resistance]

    param A,B,C,D: 4 arrays containing the corresponding values for the transfer matricies of each frequency.

    param parameter_def: the definition of the Iout parameter, indicating how it ought to be formatted

    returns: two columns (2 by n_frequencies array), and returns the parameter labels [left,right]
             and the unit labels [left,right]
    '''

    Vs = VsRsRl[0]# Thevanin source, Thevanin resistance, Load resistance
    Rs = VsRsRl[1]
    Rl = VsRsRl[2]
    n_freq = len(A)
    # Create two empty columns for real and imaginary values, will format and scale later.
    Re_Iout = np.zeros(n_freq,dtype = complex)
    Im_Iout = np.zeros(n_freq,dtype = complex)
    # Create the labels for the parameter and unit
    label_parameter = "Iout"
    label_unit      = "A"
    # For every freequency calculate the Iout value: Iin / (C*Rl + D)
    for i in range(n_freq):       
        Zin = (A[i]*Rl + B[i]) / (C[i]*Rl + D[i])
        Iin = (Vs)/(Zin + Rs)
        Iout = (Iin) / (C[i]*Rl + D[i])
        Re_Iout[i] = Iout.real
        Im_Iout[i] = Iout.imag

    # Now check for prefixes and scale
    match = re.search(r'pA', parameter_def)
    if match:
        Re_Iout = Re_Iout * 10**12
        Im_Iout = Im_Iout * 10**12
        label_unit = "p" + label_unit
    match = re.search(r'nA', parameter_def)
    if match:
        Re_Iout = Re_Iout * 10**9
        Im_Iout = Im_Iout * 10**9
        label_unit = "n" + label_unit
    match = re.search(r'uA', parameter_def)
    if match:
        Re_Iout = Re_Iout * 10**6
        Im_Iout = Im_Iout * 10**6
        label_unit = "u" + label_unit
    match = re.search(r'mA', parameter_def)
    if match:
        Re_Iout = Re_Iout * 10**3
        Im_Iout = Im_Iout * 10**3
        label_unit = "m" + label_unit
    match = re.search(r'kA', parameter_def)
    if match:
        Re_Iout = Re_Iout * 10**-3
        Im_Iout = Im_Iout * 10**-3
        label_unit = "k" + label_unit
    match = re.search(r'MA', parameter_def)
    if match:
        Re_Iout = Re_Iout * 10**-6
        Im_Iout = Im_Iout * 10**-6
        label_unit = "M" + label_unit
    match = re.search(r'GA', parameter_def)
    if match:
        Re_Iout = Re_Iout * 10**-9
        Im_Iout = Im_Iout * 10**-9
        label_unit = "G" + label_unit

    # If there is a dB prefix before the A unit, convert to the dB format and amend labels,
    match = re.search(r'dB.?A', parameter_def)
    if match:
        complex_values = Re_Iout + Im_Iout * 1j   
        out_left = 20 * np.log10( np.abs(complex_values) )       # Populate output columns correctly. (to find dB: 10log10( abs of complex value ) )
        [max(x, -160) for x in out_left]                            # Limiting dB output to be greater than -160
        out_right = np.angle(complex_values)
        label_parameter_left    = "|" + label_parameter + "|"       # Compose correct labels for the dB case
        label_parameter_right   = "/_" + label_parameter
        label_unit_left         = "dB" + label_unit
        label_unit_right        = "Rads"
    else:
        out_left = Re_Iout       # Output columns
        out_right = Im_Iout
        label_parameter_left    = "Re(" + label_parameter + ")"        # Compose correct labels for normal case
        label_parameter_right   = "Im(" + label_parameter + ")"
        label_unit_left         = label_unit
        label_unit_right        = label_unit
    
    labels = [label_parameter_left,label_parameter_right]   # Combine the labels and units into two arrays for simplicity
    units  = [label_unit_left,label_unit_right]
    out_left = out_left.real                                # Return just the real parts of the values
    out_right = out_right.real
    
    sig_figs = 4    # Formats the values to be to 4 significant figures and in exponent form.
    out_left= [f"{x:.{sig_figs-1}e}" for x in out_left]
    out_right = [f"{x:.{sig_figs-1}e}" for x in out_right]
    
    return out_left,out_right,labels,units

def calc_Pout(VsRsRl, A,B,C,D, parameter_def):
    '''
    Using the transfer function matricies of the cascaded circuit as well as source and load circuits to calculate the values
    for Pout. These are formatted and output as 2 columns.

    param VsRsRl: Array containing [Source voltage, source resistance, load resistance]

    param A,B,C,D: 4 arrays containing the corresponding values for the transfer matricies of each frequency.

    param parameter_def: the definition of the Pout parameter, indicating how it ought to be formatted

    returns: two columns (2 by n_frequencies array), and returns the parameter labels [left,right]
             and the unit labels [left,right]
    '''

    Vs = VsRsRl[0]# Thevanin source, Thevanin resistance, Load resistance
    Rs = VsRsRl[1]
    Rl = VsRsRl[2]
    n_freq = len(A)
    # Create two empty columns for real and imaginary values, will format and scale later.
    Re_Pout = np.zeros(n_freq,dtype = complex)
    Im_Pout = np.zeros(n_freq,dtype = complex)
    # Create the labels for the parameter and unit
    label_parameter = "Pout"
    label_unit      = "W"
    # For every freequency calculate the Pout value: Pin * Ap (Pin copied from above, Ap = Av* conj(Ai) )
    for i in range(n_freq):        
        Ai = (1) / (C[i]*Rl + D[i])
        Av = (Rl) / (A[i]*Rl + B[i])
        Ap = Av * np.conj(Ai)

        Zin = (A[i]*Rl + B[i]) / (C[i]*Rl + D[i])
        Iin = (Vs)/(Zin + Rs)
        Vin = (Zin*Vs)/(Zin + Rs)
        Pin = Vin * np.conj(Iin)
        Pout = Pin * Ap
    
        Re_Pout[i] = Pout.real
        Im_Pout[i] = Pout.imag

    # Now check for prefixes and scale
    match = re.search(r'pW', parameter_def)
    if match:
        Re_Pout = Re_Pout * 10**12
        Im_Pout = Im_Pout * 10**12
        label_unit = "p" + label_unit
    match = re.search(r'nW', parameter_def)
    if match:
        Re_Pout = Re_Pout * 10**9
        Im_Pout = Im_Pout * 10**9
        label_unit = "n" + label_unit
    match = re.search(r'uW', parameter_def)
    if match:
        Re_Pout = Re_Pout * 10**6
        Im_Pout = Im_Pout * 10**6
        label_unit = "u" + label_unit
    match = re.search(r'mW', parameter_def)
    if match:
        Re_Pout = Re_Pout * 10**3
        Im_Pout = Im_Pout * 10**3
        label_unit = "m" + label_unit
    match = re.search(r'kW', parameter_def)
    if match:
        Re_Pout = Re_Pout * 10**-3
        Im_Pout = Im_Pout * 10**-3
        label_unit = "k" + label_unit
    match = re.search(r'MW', parameter_def)
    if match:
        Re_Pout = Re_Pout * 10**-6
        Im_Pout = Im_Pout * 10**-6
        label_unit = "M" + label_unit
    match = re.search(r'GW', parameter_def)
    if match:
        Re_Pout = Re_Pout * 10**-9
        Im_Pout = Im_Pout * 10**-9
        label_unit = "G" + label_unit

    # If there is a dB prefix before the W unit, convert to the dB format and amend labels,
    match = re.search(r'dB.?W', parameter_def)
    if match:
        complex_values = Re_Pout + Im_Pout * 1j   
        out_left = 10 * np.log10( np.abs(complex_values) )       # Populate output columns correctly. (to find dB: 10log10( abs of complex value ) )
        [max(x, -160) for x in out_left]                            # Limiting dB output to be greater than -100
        out_right = np.angle(complex_values)
        label_parameter_left    = "|" + label_parameter + "|"       # Compose correct labels for the dB case
        label_parameter_right   = "/_" + label_parameter
        label_unit_left         = "dB" + label_unit
        label_unit_right        = "Rads"
    else:
        out_left = Re_Pout       # Output columns
        out_right = Im_Pout
        label_parameter_left    = "Re(" + label_parameter + ")"        # Compose correct labels for normal case
        label_parameter_right   = "Im(" + label_parameter + ")"
        label_unit_left         = label_unit
        label_unit_right        = label_unit
    
    labels = [label_parameter_left,label_parameter_right]   # Combine the labels and units into two arrays for simplicity
    units  = [label_unit_left,label_unit_right]
    out_left = out_left.real                                # Return just the real parts of the values
    out_right = out_right.real
    
    sig_figs = 4    # Formats the values to be to 4 significant figures and in exponent form.
    out_left= [f"{x:.{sig_figs-1}e}" for x in out_left]
    out_right = [f"{x:.{sig_figs-1}e}" for x in out_right]
    
    return out_left,out_right,labels,units

def calc_Zout(VsRsRl, A,B,C,D, parameter_def):
    '''
    Using the transfer function matricies of the cascaded circuit as well as source and load circuits to calculate the values
    for Zout. These are formatted and output as 2 columns.

    param VsRsRl: Array containing [Source voltage, source resistance, load resistance]

    param A,B,C,D: 4 arrays containing the corresponding values for the transfer matricies of each frequency.

    param parameter_def: the definition of the Zout parameter, indicating how it ought to be formatted

    returns: two columns (2 by n_frequencies array), and returns the parameter labels [left,right]
             and the unit labels [left,right]
    '''

    Vs = VsRsRl[0]# Thevanin source, Thevanin resistance, Load resistance
    Rs = VsRsRl[1]
    Rl = VsRsRl[2]
    n_freq = len(A)
    # Create two empty columns for real and imaginary values, will format and scale later.
    Re_Zout = np.zeros(n_freq,dtype = complex)
    Im_Zout = np.zeros(n_freq,dtype = complex)
    # Create the labels for the parameter and unit
    label_parameter = "Zout"
    label_unit      = "Ohms"
    # For every freequency calculate the Zout value: Vout * conj(Iout) 
    for i in range(n_freq):        
        Zout = (D[i]*Rs + B[i]) / (C[i]*Rs + A[i])
        Re_Zout[i] = Zout.real
        Im_Zout[i] = Zout.imag

    # Now check for prefixes and scale
    match = re.search(r'pOhms', parameter_def)
    if match:
        Re_Zout = Re_Zout * 10**12
        Im_Zout = Im_Zout * 10**12
        label_unit = "p" + label_unit
    match = re.search(r'nOhms', parameter_def)
    if match:
        Re_Zout = Re_Zout * 10**9
        Im_Zout = Im_Zout * 10**9
        label_unit = "n" + label_unit
    match = re.search(r'uOhms', parameter_def)
    if match:
        Re_Zout = Re_Zout * 10**6
        Im_Zout = Im_Zout * 10**6
        label_unit = "u" + label_unit
    match = re.search(r'mOhms', parameter_def)
    if match:
        Re_Zout = Re_Zout * 10**3
        Im_Zout = Im_Zout * 10**3
        label_unit = "m" + label_unit
    match = re.search(r'kOhms', parameter_def)
    if match:
        Re_Zout = Re_Zout * 10**-3
        Im_Zout = Im_Zout * 10**-3
        label_unit = "k" + label_unit
    match = re.search(r'MOhms', parameter_def)
    if match:
        Re_Zout = Re_Zout * 10**-6
        Im_Zout = Im_Zout * 10**-6
        label_unit = "M" + label_unit
    match = re.search(r'GOhms', parameter_def)
    if match:
        Re_Zout = Re_Zout * 10**-9
        Im_Zout = Im_Zout * 10**-9
        label_unit = "G" + label_unit

    # If there is a dB prefix before the Ohms unit, convert to the dB format and amend labels,
    match = re.search(r'dB.?Ohms', parameter_def)
    if match:
        complex_values = Re_Zout + Im_Zout * 1j   
        out_left = 20 * np.log10( np.abs(complex_values) )       # Populate output columns correctly. (to find dB: 10log10( abs of complex value ) )
        [max(x, -160) for x in out_left]                            # Limiting dB output to be greater than -160
        out_right = np.angle(complex_values)
        label_parameter_left    = "|" + label_parameter + "|"       # Compose correct labels for the dB case
        label_parameter_right   = "/_" + label_parameter
        label_unit_left         = "dB" + label_unit
        label_unit_right        = "Rads"
    else:
        out_left = Re_Zout       # Output columns
        out_right = Im_Zout
        label_parameter_left    = "Re(" + label_parameter + ")"        # Compose correct labels for normal case
        label_parameter_right   = "Im(" + label_parameter + ")"
        label_unit_left         = label_unit
        label_unit_right        = label_unit
    
    labels = [label_parameter_left,label_parameter_right]   # Combine the labels and units into two arrays for simplicity
    units  = [label_unit_left,label_unit_right]
    out_left = out_left.real                                # Return just the real parts of the values
    out_right = out_right.real
    
    sig_figs = 4    # Formats the values to be to 4 significant figures and in exponent form.
    out_left= [f"{x:.{sig_figs-1}e}" for x in out_left]
    out_right = [f"{x:.{sig_figs-1}e}" for x in out_right]
    
    return out_left,out_right,labels,units


def calc_Av(VsRsRl, A,B,C,D, parameter_def):
    '''
    Using the transfer function matricies of the cascaded circuit as well as source and load circuits to calculate the values
    for Av. These are formatted and output as 2 columns.

    param VsRsRl: Array containing [Source voltage, source resistance, load resistance]

    param A,B,C,D: 4 arrays containing the corresponding values for the transfer matricies of each frequency.

    param parameter_def: the definition of the Av parameter, indicating how it ought to be formatted

    returns: two columns (2 by n_frequencies array), and returns the parameter labels [left,right]
             and the unit labels [left,right]
    '''

    Vs = VsRsRl[0]# Thevanin source, Thevanin resistance, Load resistance
    Rs = VsRsRl[1]
    Rl = VsRsRl[2]
    n_freq = len(A)
    # Create two empty columns for real and imaginary values, will format and scale later.
    Re_Av = np.zeros(n_freq,dtype = complex)
    Im_Av = np.zeros(n_freq,dtype = complex)
    # Create the labels for the parameter and unit
    label_parameter = "Av"
    label_unit      = ""
    # For every freequency calculate the Av value: Vout * conj(Iout) 
    for i in range(n_freq):        
        Av = (Rl) / (A[i]*Rl + B[i])
        Re_Av[i] = Av.real
        Im_Av[i] = Av.imag


    # Bespoke way to find if there is a prefix since the unit L wont be present, but a prefix might still be requested
    if re.search(r'p', parameter_def):
        Re_Av = Re_Av * 10**12
        Im_Av = Im_Av * 10**12
        label_unit = "p" 
    elif re.search(r'n', parameter_def):
        Re_Av = Re_Av * 10**9
        Im_Av = Im_Av * 10**9
        label_unit = "n" 
    elif re.search(r'u', parameter_def):
        Re_Av = Re_Av * 10**6
        Im_Av = Im_Av * 10**6
        label_unit = "u"
    elif re.search(r'm', parameter_def):
        Re_Av = Re_Av * 10**3
        Im_Av = Im_Av * 10**3
        label_unit = "m" 
    elif re.search(r'k', parameter_def):
        Re_Av = Re_Av * 10**-3
        Im_Av = Im_Av * 10**-3
        label_unit = "k" 
    elif re.search(r'M', parameter_def):
        Re_Av = Re_Av * 10**-6
        Im_Av = Im_Av * 10**-6
        label_unit = "M"
    elif re.search(r'G', parameter_def):
        Re_Av = Re_Av * 10**-9
        Im_Av = Im_Av * 10**-9
        label_unit = "G"

    # If there is a dB prefix before the L unit, convert to the dB format and amend labels,
    match = re.search(r'dB', parameter_def)
    if match:
        complex_values = Re_Av + Im_Av * 1j   
        out_left = 20 * np.log10( np.abs(complex_values) )       # Populate output columns correctly. (to find dB: 10log10( abs of complex value ) )
        [max(x, -160) for x in out_left]                            # Limiting dB output to be greater than -160
        out_right = np.angle(complex_values)
        label_parameter_left    = "|" + label_parameter + "|"       # Compose correct labels for the dB case
        label_parameter_right   = "/_" + label_parameter
        label_unit_left         = "dB" + label_unit
        label_unit_right        = "Rads"
    else:
        out_left = Re_Av       # Output columns
        out_right = Im_Av
        label_parameter_left    = "Re(" + label_parameter + ")"        # Compose correct labels for normal case
        label_parameter_right   = "Im(" + label_parameter + ")"
        label_unit_left         = label_unit + "L"
        label_unit_right        = label_unit + "L"
    
    labels = [label_parameter_left,label_parameter_right]   # Combine the labels and units into two arrays for simplicity
    units  = [label_unit_left,label_unit_right]
    out_left = out_left.real                                # Return just the real parts of the values
    out_right = out_right.real
    
    sig_figs = 4    # Formats the values to be to 4 significant figures and in exponent form.
    out_left= [f"{x:.{sig_figs-1}e}" for x in out_left]
    out_right = [f"{x:.{sig_figs-1}e}" for x in out_right]
    

    return out_left,out_right,labels,units


def calc_Ai(VsRsRl, A,B,C,D, parameter_def):
    '''
    Using the transfer function matricies of the cascaded circuit as well as source and load circuits to calculate the values
    for Ai. These are formatted and output as 2 columns.

    param VsRsRl: Array containing [Source voltage, source resistance, load resistance]

    param A,B,C,D: 4 arrays containing the corresponding values for the transfer matricies of each frequency.

    param parameter_def: the definition of the Ai parameter, indicating how it ought to be formatted

    returns: two columns (2 by n_frequencies array), and returns the parameter labels [left,right]
             and the unit labels [left,right]
    '''

    Vs = VsRsRl[0]# Thevanin source, Thevanin resistance, Load resistance
    Rs = VsRsRl[1]
    Rl = VsRsRl[2]
    n_freq = len(A)
    # Create two empty columns for real and imaginary values, will format and scale later.
    Re_Ai = np.zeros(n_freq,dtype = complex)
    Im_Ai = np.zeros(n_freq,dtype = complex)
    # Create the labels for the parameter and unit
    label_parameter = "Ai"
    label_unit      = ""
    # For every freequency calculate the Ai value: Vout * conj(Iout) 
    for i in range(n_freq):        
        Ai = (1) / (C[i]*Rl + D[i])
        Re_Ai[i] = Ai.real
        Im_Ai[i] = Ai.imag

    # Bespoke way to find if there is a prefix since the unit L wont be present, but a prefix might still be requested
    if re.search(r'p', parameter_def):
        Re_Ai = Re_Ai * 10**12
        Im_Ai = Im_Ai * 10**12
        label_unit = "p" 
    elif re.search(r'n', parameter_def):
        Re_Ai = Re_Ai * 10**9
        Im_Ai = Im_Ai * 10**9
        label_unit = "n" 
    elif re.search(r'u', parameter_def):
        Re_Ai = Re_Ai * 10**6
        Im_Ai = Im_Ai * 10**6
        label_unit = "u"
    elif re.search(r'm', parameter_def):
        Re_Ai = Re_Ai * 10**3
        Im_Ai = Im_Ai * 10**3
        label_unit = "m" 
    elif re.search(r'k', parameter_def):
        Re_Ai = Re_Ai * 10**-3
        Im_Ai = Im_Ai * 10**-3
        label_unit = "k" 
    elif re.search(r'M', parameter_def):
        Re_Ai = Re_Ai * 10**-6
        Im_Ai = Im_Ai * 10**-6
        label_unit = "M"
    elif re.search(r'G', parameter_def):
        Re_Ai = Re_Ai * 10**-9
        Im_Ai = Im_Ai * 10**-9
        label_unit = "G"

    # If there is a dB prefix, convert to the dB format and amend labels,
    match = re.search(r'dB', parameter_def)
    if match:
        complex_values = Re_Ai + Im_Ai * 1j   
        out_left = 20 * np.log10( np.abs(complex_values) )       # Populate output columns correctly. (to find dB: 10log10( abs of complex value ) )
        [max(x, -160) for x in out_left]                            # Limiting dB output to be greater than -160
        out_right = np.angle(complex_values)
        label_parameter_left    = "|" + label_parameter + "|"       # Compose correct labels for the dB case
        label_parameter_right   = "/_" + label_parameter
        label_unit_left         = "dB" + label_unit
        label_unit_right        = "Rads"
    else:
        out_left = Re_Ai       # Output columns
        out_right = Im_Ai
        label_parameter_left    = "Re(" + label_parameter + ")"        # Compose correct labels for normal case
        label_parameter_right   = "Im(" + label_parameter + ")"
        label_unit_left         = label_unit + "L"
        label_unit_right        = label_unit + "L"
    
    labels = [label_parameter_left,label_parameter_right]   # Combine the labels and units into two arrays for simplicity
    units  = [label_unit_left,label_unit_right]
    out_left = out_left.real                                # Return just the real parts of the values
    out_right = out_right.real
        
    sig_figs = 4    # Formats the values to be to 4 significant figures and in exponent form.
    out_left= [f"{x:.{sig_figs-1}e}" for x in out_left]
    out_right = [f"{x:.{sig_figs-1}e}" for x in out_right]
    
    return out_left,out_right,labels,units

#===========================================================================================================================================================================================================
#================ MAIN SCRIPT ==============================================================================================================================================================================
#===========================================================================================================================================================================================================

                                    # Defining command line arguments 
parser = argparse.ArgumentParser()
# First argument is input file:
parser.add_argument('input_file', metavar='INPUT_FILE', nargs='?', type=str, help='input .net file')   
# Second argument is output file:
parser.add_argument('output_file', metavar='OUTPUT_FILE', nargs='?', type=str, help='output .csv file')
# -i: Takes in string
parser.add_argument('-i', metavar='NAME', type=str, help='name of input and output files')              
# -p: Takes in multiple ints
parser.add_argument('-p', metavar='PARAMS', type=int, nargs='+', help='columns to plot')                
args = parser.parse_args()


# Determine input and output file names (if -i provided use that, else use the ones provided)
if args.i is not None:
    input_file = args.i + '.net'
    output_file = args.i + '.csv'
else:
    input_file = args.input_file
    output_file = args.output_file

# Check if input and output file names were provided
if input_file is None or output_file is None:
    print('Please provide input and output file names.')
    sys.exit()


with open(output_file, 'w') as creating_new_csv_file: 
   pass # Create empty file.

# Open and read .net file
net_file    = open(input_file,'rt')
read_file   = net_file.readlines()

# Calling Functional Blocks
circuit_cascaded            = cascaded_circuit(read_file)

frequency_points            = frequencies(read_file)

VsRsRl                      = get_VsRsRl(read_file)

[A,B,C,D]                   = transfer_functions(circuit_cascaded,frequency_points)

built_array,parameters      = output_array(read_file, frequency_points)

finalised_array = calculations(A,B,C,D, built_array, parameters, VsRsRl)

# This adds an extra comma to the end of the value lines, to exactly conform with model outputs.
for i in range(len(frequency_points)):
    finalised_array[i+2].append("") 

# Write the finalised array to a .csv file
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(finalised_array)


# Use array of numbers if provided
if args.p is not None:
    print_columns = args.p
    # Plot the columns requested if -p argument present
    for i in range(len(print_columns)):

        column = print_columns[i]           # If a column number which is too large is provided, return an error
        if i > len(finalised_array[0]):
            raise Exception("Column number provided is out of bounds.")
    
        x = []  # Fill out x and y points by looking at the finalised output array
        y = []
        for j in range(len(finalised_array) - 2):
            x.append(finalised_array[2+j][0])
            y.append(finalised_array[2+j][column])

        type = finalised_array[0][1+i]  # Retrieve the plotted value type and unit
        unit = finalised_array[1][1+i]

        plt.figure()    
        plt.plot(x,y,'-o')      # Plot the points, marks are circles
        plt.yticks(rotation=45)     # Rotate y tick labels, more space for the axis label 
        plt.xlabel('Frequency (Hz)')    # x-axis label
        plt.ylabel( type.strip() + ' (' + unit.strip() + ')')   # y-axis label, take from the top two rows and strip
        plt.tick_params(axis='both', which='major', labelsize=6)    # Make tick labels smaller for more space and less overlay
         # Save it to a png, same name as the output file with column number appended, increase resolution
        plt.savefig(output_file[:-4]+ "_" + str(print_columns[i]) + ".png", dpi=600)       

