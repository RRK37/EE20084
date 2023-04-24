####################################################################################################################################################################
#       Filename:       mainScript.py
#       Summary:        This script simulates a circuit provided as a .net file, outputting the results of the simulation in a .csv file, and possibly graph plots.
#       Description:    Taking in a .csv file it is read and a corresponding virtual cascaded circuit is composed. Its charecteristic two-port ABCD matrix is calc- 
#                       ulated for each frequency defined in the .net file. It is then simulated in tangent to the input and output circuits provided in the .net.
#                       The results of the simulation are stored in a matrix and output in a .csv.
#       Author:         R. Klotins
####################################################################################################################################################################


#=========== Libraries  ============================================================================================================================================
import re



#=========== Accesiblity functions =================================================================================================================================

def get_n1(component_def):
    '''
    Returns the n1 integer value of a string component definition. If it cannot be found it returns -1.

    param component_def: the component definition (eg. "n1=1 n2=2 R=3")
    
    returns: integer value of n1 (eg. 1)    
    '''
    match = re.search(r'n1\s*=\s*(\d+)', component_def)
    if match:
        return int(match.group(1))
    else:
        return -1


def get_n2(component_def):
    '''
    Returns the n2 integer value of a string component definition. If it cannot be found it returns -1.

    param component_def: the component definition (eg. "n1=1 n2=2 R=3")
    
    returns: integer value of n2 (eg. 2)    
    '''
   
    match = re.search(r'n2\s*=\s*(\d+)', component_def)
    if match:
        return int(match.group(1))
    else:
        return -1
    
def get_val(component_def):
    '''
    Returns the value and type of the component definition provided. The type can either be
    1: R, 2: G, 3: L, 4: C. If not found, returns -1.

    param component_def: the component definition (eg. "n1=1 n2=2 C=3")
    
    returns: integer value of n1 in the first position of vector, and the type in the second (eg. [3,4])     
    '''

    component_type = -1 # If none of the definition types are matched, return -1 (the "not matched flag")
    component_val = -1

    match = re.search(r"R\s*=\s*(\d+)", component_def)  # Resistance: component type 1
    if match:
        component_val = int(match.group(1))
        component_type = 1
    match = re.search(r"G\s*=\s*(\d+)", component_def)  # Conductance: component type 2
    if match:
        component_val = int(match.group(1))
        component_type = 2
    match = re.search(r"L\s*=\s*(\d+)", component_def)  # Inductance: component type 3
    if match:
        component_val = int(match.group(1))
        component_type = 3
    match = re.search(r"C\s*=\s*(\d+)", component_def)  # Capacitance: component type 4
    if match:
        component_val = int(match.group(1))
        component_type = 4

    return [component_val,component_type]   


#=========== Functional Block 1: Cascaded Circuit Composition ======================================================================================================

def cascaded_circuit(file):
    '''
    Takes in the .net file in string array from, every line being in a new index. Retruns the 
    component defintions in a string array, sorted according to how they would appear in a 
    cascaded circuit diagram. 
    
    
    
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

