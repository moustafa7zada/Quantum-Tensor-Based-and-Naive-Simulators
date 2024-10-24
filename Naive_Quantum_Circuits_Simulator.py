import numpy as np 

class Naive_Circuit() : 
    def __init__(self , num_of_qubits:int , keep_already_made_gates = False ) : 
        if num_of_qubits == 0 : 
            raise Exception("You Cannot Have a Circuit wiht No qubits at ALL!")
        self.keep_already_made_gates = keep_already_made_gates
        
        if keep_already_made_gates : 
            self.already_made_gates = {}
            
        self.num_of_qubits = num_of_qubits 
        self.statevec = np.zeros(2**num_of_qubits)
        self.statevec[0] = 1 
        self.sequence = []    
        self.identity = np.array(([1,0],[0,1]))
        self.X_matrix = np.array(([0,1],[1,0]))
        self.Z_matrix = np.array(([1,0], [0,-1] ))
        self.H_matrix = np.array(([1.0/np.sqrt(2) , 1.0/np.sqrt(2)] , [1.0/np.sqrt(2) , -1.0/np.sqrt(2)]))
        self.CNOT_matrix = np.array(([1,0,0,0],
                                     [0,1,0,0],
                                     [0,0,0,1],
                                     [0,0,1,0]      ))
        
        self.zero_proj = np.array(([1,0] ,
                                   [0,0]))
        self.one_proj = np.array(([0,0], 
                                  [0,1]))
        
    def Measure(self) : 
        return f'{np.random.choice(range(len(self.statevec)),p=np.square(self.statevec)):#0{self.num_of_qubits+2}b}'[2:]
    
    
    def ExpVal(self) : 
        big_Z = self.Z_matrix
        for i in range(self.num_of_qubits - 1) : 
            big_Z = np.kron(big_Z  , self.Z_matrix)
        return np.dot(np.dot( self.statevec , big_Z ) ,self.statevec.T.conjugate() )
    
    def reset(self) : 
        self.statevec = np.zeros(2**self.num_of_qubits)
        self.statevec[0] = 1
        self.sequence = []        
        print("Your Circuit Has Been Reset!")
        
    def draw(self) : 
        print("A Simple Representation of Your Circuit \n \t" ,  self.sequence)
        
        
    def X(self , position:int) : 
        if position >= self.num_of_qubits or position < 0 : 
            raise Exception(f"Position of the X gate is not Correct, you only have {self.num_of_qubits} qubits and you tried applying it to the {position+1}th qubit")
            
        #lets check if we have already made such a gate , if not, lets make an X gate that can be mulitplied by our whole StateVector
        try : 
            self.statevec = np.dot(self.already_made_gates[f'X on the {position}'] , self.statevec) 
            
        except : 
            The_Scaled_UP_Matrix = None 
            for i in range(self.num_of_qubits) : 
                if The_Scaled_UP_Matrix is None and i != position : 
                    The_Scaled_UP_Matrix = self.identity
                    continue
                if The_Scaled_UP_Matrix is None :
                    The_Scaled_UP_Matrix = self.X_matrix
                    continue
                if i == position : 
                    The_Scaled_UP_Matrix = np.kron(The_Scaled_UP_Matrix , self.X_matrix ) 
                else : 
                    The_Scaled_UP_Matrix = np.kron(The_Scaled_UP_Matrix , self.identity )
                
            self.statevec = np.dot(The_Scaled_UP_Matrix , self.statevec)  
            if self.keep_already_made_gates : 
                self.already_made_gates[f'X on the {position}'] = The_Scaled_UP_Matrix
        self.sequence.append(f"X[{position}]")
            
                
    def Z(self , position:int) : 
        if position >= self.num_of_qubits or position < 0 : 
            raise Exception(f"Position of the Z gate is not correct, you only have {self.num_of_qubits} qubits and you tried applying it to the {position+1}th qubit")
            
        #lets make a Z gate that can be mulitplied by our whole StateVector
        try : 
            self.statevec= np.dot(self.already_made_gates[f'Z on the {position}'] , self.statevec) 
        except : 
            The_Scaled_UP_Matrix = None 
            for i in range(self.num_of_qubits) : 
                if The_Scaled_UP_Matrix is None and i != position : 
                    The_Scaled_UP_Matrix = self.identity
                    continue
                if The_Scaled_UP_Matrix is None :
                    The_Scaled_UP_Matrix = self.Z_matrix
                    continue
                if i == position : 
                    The_Scaled_UP_Matrix = np.kron(The_Scaled_UP_Matrix , self.Z_matrix ) 
                else : 
                    The_Scaled_UP_Matrix = np.kron(The_Scaled_UP_Matrix , self.identity )
                
                
            self.statevec = np.dot(The_Scaled_UP_Matrix , self.statevec)
            if self.keep_already_made_gates : 
                 self.already_made_gates[f'Z on the {position}'] = The_Scaled_UP_Matrix

        self.sequence.append(f"Z[{position}]")


    def H(self , position:int) : 
        if position >= self.num_of_qubits or position < 0 : 
            raise Exception(f"Position of the Hadamard gate is not correct, you only have {self.num_of_qubits} qubits and you tried applying it to the {position+1}th qubit")
            
        #lets make an H gate that can be mulitplied by our whole StateVector
        try : 
            self.statevec= np.dot(self.already_made_gates[f'H on the {position}'] , self.statevec) 

        except:
            The_Scaled_UP_Matrix = None 
            for i in range(self.num_of_qubits) : 
                if The_Scaled_UP_Matrix is None and i != position : 
                    The_Scaled_UP_Matrix = self.identity
                    continue
                if The_Scaled_UP_Matrix is None :
                    The_Scaled_UP_Matrix = self.H_matrix
                    continue
                if i == position : 
                    The_Scaled_UP_Matrix = np.kron(The_Scaled_UP_Matrix , self.H_matrix ) 
                else : 
                    The_Scaled_UP_Matrix = np.kron(The_Scaled_UP_Matrix , self.identity )
                
            self.statevec = np.dot(The_Scaled_UP_Matrix , self.statevec)
            if self.keep_already_made_gates : 
                self.already_made_gates[f'H on the {position}'] = The_Scaled_UP_Matrix

        self.sequence.append(f"Hadamard[{position}]")
        
    def I(self , position:int) : 
        if position >= self.num_of_qubits or position < 0  : 
            raise Exception(f"Position of the identity gate is not correct, you only have {self.num_of_qubits} qubits and you tried applying it to the {position+1}th qubit")
        self.sequence.append(f"X[{position}]")

            
            
            
    def CNOT(self , control:int ,target:int ) : 
        if control >= self.num_of_qubits or target >= self.num_of_qubits or control < 0 or target < 0 : 
            raise Exception(f"Position of the CNOT gate is not correct, you only have {self.num_of_qubits} qubits and you tried applying it to a qubit that you dont have either as the target or as the control qubit ")
        
        try : 
            self.statevec= np.dot(self.already_made_gates[f'CNOT from {control} to {target}'] , self.statevec) 

        except : 
            #lets make a CNOT gate that can be mulitplied by our whole StateVector
            The_Scaled_UP_Matrix_zero = None 
            The_Scaled_UP_Matrix_one = None 
            for i in range(self.num_of_qubits) : 
                if The_Scaled_UP_Matrix_one is None and i == control : 
                    The_Scaled_UP_Matrix_zero = self.zero_proj
                    The_Scaled_UP_Matrix_one = self.one_proj
                    continue
                if The_Scaled_UP_Matrix_one is None and i == target :  
                    The_Scaled_UP_Matrix_zero = self.identity
                    The_Scaled_UP_Matrix_one = self.X_matrix
                    continue
                if The_Scaled_UP_Matrix_one is None : 
                    The_Scaled_UP_Matrix_zero = self.identity
                    The_Scaled_UP_Matrix_one = self.identity
                    continue
                if i != target : 
                    if i == control : 
                        The_Scaled_UP_Matrix_zero = np.kron(The_Scaled_UP_Matrix_zero , self.zero_proj) 
                        The_Scaled_UP_Matrix_one = np.kron(The_Scaled_UP_Matrix_one , self.one_proj)
                        continue
                    The_Scaled_UP_Matrix_zero = np.kron(The_Scaled_UP_Matrix_zero , self.identity)
                    The_Scaled_UP_Matrix_one = np.kron(The_Scaled_UP_Matrix_one , self.identity ) 
                    continue
                else : 
                    The_Scaled_UP_Matrix_zero = np.kron(The_Scaled_UP_Matrix_zero , self.identity)
                    The_Scaled_UP_Matrix_one = np.kron(The_Scaled_UP_Matrix_one , self.X_matrix )
                    continue
        
            The_Scaled_UP_Matrix_zero += The_Scaled_UP_Matrix_one
            self.statevec = np.dot(The_Scaled_UP_Matrix_zero , self.statevec)  
            if self.keep_already_made_gates : 
                self.already_made_gates[f'CNOT from {control} to {target}'] = The_Scaled_UP_Matrix_zero
            
        self.sequence.append(f"CNOT[{control} --> {target}")