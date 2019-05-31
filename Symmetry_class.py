import re as re
import numpy as np

class Symmetry:
    def __init__(self, matrix, operation):
        self.matrix = matrix
        self.operation = operation
        self.char = np.trace(matrix)
    def __str__(self):
        return '{} {}\n'.format(self.matrix,self.operation)
    def __repr__(self):
        return '{} {}\n'.format(self.matrix,self.operation)

class SymmClass:
    def __init__(self, charachter, operation):
        self.charachter =charachter
        self.operation =operation
    def __str__(self):
        return '{} {}\n'.format(self.charachter,self.operation)
    def __repr__(self):
        return '{} {}\n'.format(self.charachter,self.operation) 