"""General Huckel Solver.

For an arbritary linear, cyclic or a sp2 Platonic polyene consisting of a 
carbon skeleton with n atoms, this script will generate for the π-system the 
Huckel energies and their associated degeneracies. The octahedral shape is also
supported. It also runs a basic check that the energies are correct by summing 
the Huckel energies together, which should equal 0. The energies and 
degeneracies are directly printed as well as saved in a .csv file.

This script requires the import of the "numpy", "decimal", "csv" and "sys"
modules.

The program will require user input at various stages. 
To determine the kind of molecule you wish to calculate the energies for, 
please type "linear" for a linear molecule, "cyclic" for a cyclic molecule or 
"Platonic solid" for a sp2 hybridized Platonic solid (capitalization does 
not matter) when prompted.

To determine the number of carbons in your linear or cyclic system, please 
type an integer number when prompted. 

To choose the kind of Platonic solid you wish to solve, please type 
"tetrahedron" for a tetrahedron, "cube" for a cube, "octahedron" for an 
octahedron and "dodecahedron" for a dodecahedron when prompted.
"""

import sys
import numpy as np
from decimal import Decimal
import csv

matrix = []
alpha = 0
beta = -1

def General_Huckel_Solver():
    """Order the script functions to form a coherent Huckel solver script."""
    
    global matrix
    global alpha
    global beta
    
    while True:             
        problem_type = input("Do you have a linear polyene, cyclic polyene or" 
                    " a sp2 hybridized Platonic solid? ")
        if problem_type.lower() == "linear":
            matrix = Huckel_Solver_linear_n_polyene()
            break
        elif problem_type.lower() == "cyclic":
            matrix = Huckel_Solver_cyclic_n_polyene()
            break
        elif problem_type.lower() == "platonic solid":
            matrix = Huckel_Solver_Platonic_Solid()
            break
        elif problem_type == "exit":
            sys.exit()
        else:
            print()
            print('You should input "linear", "cyclic" or "Platonic solid." '
                  'To exit the program, type exit.')
    
    eval_list_fun = get_evals()
    
    # Determines the sum of the eigenvalues - this should be 0 under Huckel 
    # approximations.
    eigenval_sum = "%.3f" % sum(eval_list_fun)
    evalue_sum_statement = "The sum of the Huckel energies is " + eigenval_sum
    
    final_list = Degeneracy_count()
    
    print()
    for i in final_list:
        print(i)
    print()
    print(evalue_sum_statement)
    
    # Writes the Huckel energy levels and their degeneracies to a .csv file.
    f = open("Huckel_Energies.csv", 'w', newline='')
    for item in Degeneracy_count():
        csv.writer(f).writerow(item)
    csv.writer(f).writerow([evalue_sum_statement])   
    f.close()
    


def Huckel_Solver_linear_n_polyene():
    """Form the Huckel matrix for a linear polyene with n carbons."""

    # Exception handling to ensure that when a physically insignificant number
    # of carbons (non-integer) are inputted, the program does not crash.
    while True:
        print()
        num = input("How many carbons in your linear poly-ene? ")
        try:
            n = int(num)
            break
        except ValueError:
            print()
            print("Please enter an integer number of carbons.")
        
    base_matrix_linear = np.zeros((n, n))
    
    # Add in the alpha-interactions for the leading diagonal.
    for i in range(n):
        base_matrix_linear[i,i] = alpha
    
    # Add in the beta-interactions for the subdiagonal.
    for i in range(1,n,1):
        base_matrix_linear[i,i-1] = beta 
    
    # Add in the beta-interactions for the superdiagonal.
    for i in range(1, n, 1):
        base_matrix_linear[i-1,i] = beta 
        
    return base_matrix_linear



def Huckel_Solver_cyclic_n_polyene():
    """Form the Huckel matrix for a cyclic polyene with n carbons."""
    
    # Exception handling to ensure that when a physically insignificant number
    # of carbons (non-integer) are inputted, the program does not crash.
    while True:
        print()
        num = input("How many carbons in your cyclic poly-ene? ")
        try:
            n = int(num)
            break
        except ValueError:
            print()
            print("Please enter an integer number of carbons.")
            
    base_matrix_cyclic = np.zeros((n, n))
    
    # Add in the alpha-interactions for the leading diagonal.
    for i in range(n):
        base_matrix_cyclic[i,i] = alpha
    
    # Add in the beta-interactions for the subdiagonal.
    for i in range(1,n,1):
        base_matrix_cyclic[i,i-1] = beta 
    
    # Add in the beta-interactions for the superdiagonal.
    for i in range(1, n, 1):
        base_matrix_cyclic[i-1,i] = beta 
        
    # Add in the first (1-N) cyclic beta interaction.
    base_matrix_cyclic[0,n-1] = beta

    # Add in the second (N-1) cyclic beta interaction.
    base_matrix_cyclic[n-1,0] = beta

    return base_matrix_cyclic




def get_evals():
    """ Return the eigenvalues of a matrix."""
    
    evalues= np.linalg.eigvals(matrix)
    
    # Sorts the eigenvalues from lowest to highest.
    evalues_numpy = np.sort(evalues)
    evalues_list = list(evalues_numpy)
    
    return evalues_list


def Degeneracy_count():
    """ Work out the degeneracy of each eigenvalue in a given list."""
    eval_list = get_evals()
    
    eigenval_sum = "%.3f" % sum(eval_list)
    evalue_sum_statement = "The sum of the Huckel energies is "
    evalue_sum_list = [evalue_sum_statement, eigenval_sum]

   
    # Rounds the eigenvalues to 3 significant figures and makes a list of 
    # eigenvalue/degeneracy pairings.
    nice_eval_list = [float(Decimal("%.3f" % elem)) for elem in eval_list]
    degeneracy_base = [["Energy", "Degeneracy"]]
    degeneracies = [[x,nice_eval_list.count(x)] for x in set(nice_eval_list)]
    degeneracies_sorted = sorted(degeneracies)
    degen_final = degeneracy_base + degeneracies_sorted
    return degen_final
    
def Huckel_Solver_Platonic_Solid():
    
    alpha = 0
    beta = -1
    
    while True:
        print()
        shape = input("What kind of shape do you have? ")
        
        # Forms the Huckel matrix for the tetrahedron.
        if shape.lower() == "tetrahedron":
            # For tetrahedron, all non-diagonal elements are the beta
            # interactions. The Huckel matrix is therefore constructed from a 
            # base matrix where all terms are beta-interactions.
            base_matrix_tetrahedron = np.array([[beta, beta, beta, beta], \
                                                [beta, beta, beta, beta], \
                                                [beta, beta, beta, beta], \
                                                [beta, beta, beta, beta]])
            # Add in the alpha-interactions for the leading diagonal.
            for i in range(4):
                base_matrix_tetrahedron[i,i] = alpha
                
            return base_matrix_tetrahedron
            
        
        # Forms the Huckel matrix for the octahedron.   
        elif shape.lower() == "octahedron":
            base_matrix_octahedron = np.zeros((6,6))
            # This changes all matrix elements to a beta-interaction. 
            # Editing of the matrix further down will remove beta-interactions
            # which are not physically present.
            huckel_matrix_octahedron = base_matrix_octahedron + beta
            
            
            for i in range(6):
                # Add in the alpha-interactions for the leading diagonal.
                huckel_matrix_octahedron[i,i] = alpha
                # These are the matrix elements from non-adjacent vertices.
                huckel_matrix_octahedron[i,6-i-1] = 0
                
            return huckel_matrix_octahedron
         
            
        # Forms the Huckel matrix for the cube.
        elif shape.lower() == "cube":
            base_matrix_cube = np.zeros((8,8))
            # Add in the alpha-interactions for the leading diagonal.
            for i in range(8):
                base_matrix_cube[i,i] = alpha
            # Add in the beta-interactions .
            column_elem = [[1,2,4], [0,3,5], [0,3,6], [1,2,7], [0,5,6], \
                           [1,4,7], [2,4,7], [3,5,6]]
            for i in range(len(column_elem)):
                base_matrix_cube[i,column_elem[i][0]] = beta
                base_matrix_cube[i,column_elem[i][1]] = beta
                base_matrix_cube[i,column_elem[i][2]] = beta  
                
            return base_matrix_cube
        
        
        # Forms the Huckel matrix for the dodecahedron.
        elif shape.lower() == "dodecahedron":
            base_matrix_dodec = np.zeros((20,20))
        # Add in the alpha-interactions for the leading diagonal.
            for i in range(20):
                base_matrix_dodec[i,i] = alpha
        #Add in the beta-interactions.
            column_elem_dodec = [[13,14,15],[4,5,12],[6,13,18],[7,14,19], \
                                 [1,10,18],[1,11,19],[2,10,15],[3,11,15], \
                                 [9,13,16],[8,14,17],[4,6,11],[5,7,10], \
                                 [1,16,17],[0,2,8],[0,3,9],[0,6,7], \
                                 [8,12,18], [9,12,19], [2,4,16], [3,5,17]]
            for i in range(len(column_elem_dodec)):
                base_matrix_dodec[i, column_elem_dodec[i][0]] = beta
                base_matrix_dodec[i, column_elem_dodec[i][1]] = beta
                base_matrix_dodec[i, column_elem_dodec[i][2]] = beta
            
            return base_matrix_dodec

        
        elif shape == "none":
            sys.exit()
            
        else:
            print("Your shape should be a tetrahedron, cube, octahedron "
                  "or dodecahedron. Type none to quit.")
        


if __name__ == '__main__':
    General_Huckel_Solver()
    
        
        
