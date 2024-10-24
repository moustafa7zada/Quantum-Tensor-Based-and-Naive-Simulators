import numpy as np 

class Tensor_Circuit() : 
    def __init__(self , num_of_qubits:int ) : 
        if num_of_qubits == 0 : 
            raise Exception("You Cannot Have a Circuit wiht No qubits at ALL!")
            
        self.num_of_qubits = num_of_qubits 
        self.statevec = np.zeros(tuple( 2 for _ in range(num_of_qubits) ) )
        self.statevec[tuple([0]*num_of_qubits)] = 1 
        self.sequence = []    
        self.identity = np.array(([1,0],[0,1]))
        self.X_matrix = np.array(([0,1],[1,0]))
        self.Z_matrix = np.array(([1,0], [0,-1] ))
        self.H_matrix = np.array(([1.0/np.sqrt(2) , 1.0/np.sqrt(2)] , [1.0/np.sqrt(2) , -1.0/np.sqrt(2)]))
        self.CNOT_matrix = np.array(([1,0,0,0],
                                     [0,1,0,0],
                                     [0,0,0,1],
                                     [0,0,1,0]      ))
        self.CNOT_matrix = np.reshape(self.CNOT_matrix , (2 ,2  ,2 , 2 ) ) 
        self.zero_proj = np.array(([1,0] ,
                                   [0,0]))
        self.one_proj = np.array(([0,0], 
                                  [0,1]))
        
    def Measure(self) : 
        return f'{np.random.choice(range(len(self.statevec.reshape(-1))),p=np.square(self.statevec.reshape(-1))):#0{self.num_of_qubits+2}b}'[2:]
    
    
    def ExpVal(self) : 
        big_Z = self.Z_matrix
        for i in range(self.num_of_qubits - 1) : 
            big_Z = np.kron(big_Z  , self.Z_matrix)
        return np.dot(np.dot( self.statevec.reshape(-1) , big_Z ) ,self.statevec.reshape(-1).T.conjugate() )
    
    def reset(self) : 
        self.statevec = np.zeros(( 2 for _ in range(self.num_of_qubits) ) )
        self.statevec[tuple([0]*self.num_of_qubits)] = 1 
        self.sequence = []        
        print("Your Circuit Has Been Reset!")
        
    def draw(self) : 
        print("A Simple Representation of Your Circuit \n \t" ,  self.sequence)
        
        
    def X(self , position:int) : 
        if position >= self.num_of_qubits or position < 0 : 
            raise Exception(f"Position of the X gate is not Correct, you only have {self.num_of_qubits} qubits and you tried applying it to the {position+1}th qubit")
        
        self.statevec = np.tensordot(self.X_matrix , self.statevec , axes=([1] , [position] ) ) 
        self.sequence.append(f"X[{position}]")
            
                
    def Z(self , position:int) : 
        if position >= self.num_of_qubits or position < 0 : 
            raise Exception(f"Position of the Z gate is not Correct, you only have {self.num_of_qubits} qubits and you tried applying it to the {position+1}th qubit")
        
        self.statevec = np.tensordot(self.Z_matrix , self.statevec , axes=([1] , [position] ) ) 
        self.sequence.append(f"Z[{position}]")


    def H(self , position:int) : 
        if position >= self.num_of_qubits or position < 0 : 
            raise Exception(f"Position of the H gate is not Correct, you only have {self.num_of_qubits} qubits and you tried applying it to the {position+1}th qubit")
        
        self.statevec = np.tensordot(self.H_matrix , self.statevec , axes=([1] , [position] ) ) 
        self.sequence.append(f"Hadamnard[{position}]")
        
    def I(self , position:int) : 
        if position >= self.num_of_qubits or position < 0  : 
            raise Exception(f"Position of the identity gate is not correct, you only have {self.num_of_qubits} qubits and you tried applying it to the {position+1}th qubit")
        self.sequence.append(f"I[{position}]")
            
            
    def CNOT(self , control:int ,target:int ) : 
        if control >= self.num_of_qubits or target >= self.num_of_qubits or control < 0 or target < 0 : 
            raise Exception(f"Position of the CNOT gate is not correct, you only have {self.num_of_qubits} qubits and you tried applying it to a qubit that you dont have either as the target or as the control qubit ")
        
        self.statevec = np.tensordot(self.CNOT_matrix , self.statevec , axes=( [2,3] , [control , target] ) ) 
        self.sequence.append(f"CNOT[{control} --> {target}")