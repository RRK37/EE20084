####################################################################################################
#   Filename:       mainScript.py
#   Summary:        Simulates a circuit from a given .net definition of one, outputting CSV data 


####################################################################################################

# Importing Libraries
import re

# ================== Functions =====================================================================

def get_n1 (component_def: str) -> int:    
    raise Exception("This function is not defined yet.")

def get_n2 (component_def: str) -> int:
    raise Exception("This function is not defined yet.")

def get_val (component_def: str) -> list[int]: 
    raise Exception("This function is not defined yet.")

def cascaded_circuit (file: str, n_lines: int) -> list[str]:
    """
    param file: string type of the input .net file 
    param n_lines: this is the number of lines in the input net file 
    return: this is a string array of component definitions, sorted as a cascaded curcuit 

    """

    # Finding where the <CIRCUIT> block starts by going line by line checking for <CIRCUIT>
    line_index = 0
    cir_index = 0
    for i in range(n_lines):
        if "<CIRCUIT>" in file[line_index:line_index+10]:
            cir_start = line_index + 10 

        line_index = file.find("\n",line_index+3)   # Find the next line start

    # Finding the end of the <CIRCUIT> block. Going line by line till: </CIRCUIT>
    for i in 


    # Populating the 'serial' and 'parallel' arrays with corresponding component defintions
    serial_comp = []
    parallel_comp = []

    more_nodes = True
    while more_nodes:
        match = re.search(r'n1\s*=\s*(\d+)', file[cir_index])
        if match =




# ================== Main Code =====================================================================

# Open and Read .net file
net_file    = open('D:\\Code\\EE20084\\EE20084\\input_files\\a_Test_Circuit_1.net','rt')
file        = net_file.read()
n_lines     = file.count("\n") 



cascaded_circuit(file,n_lines)




net_file.close()