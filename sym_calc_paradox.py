#!/usr/bin/python

from post_procesing_tools import *

path = "/home/mjocic/data_paradox/nbnd_120"
k_G=1
E=Gather_E(path,units='Hartree')
E_G = E[k_G-1]
k_points= Gather_K(path)
k_Gpoint = k_points[k_G-1]
G_G = Gather_G(path,k_G)
C_G = Gather_C(path,k_G)

symmetry_list = Gather_Symmetries(path)

uk_start=2
uk_stop=4
charachters=[]
for s in range(len(symmetry_list)):
    sym_name=symmetry_list[s].operation
    sym_matrix=symmetry_list[s].matrix
    print('\nNow performing for '+ sym_name)
    char=0
    for k in range(uk_start,uk_stop+1):
        RC_G = np.zeros(len(G_G),dtype=complex)
        for i in range(len(G_G)):
            g_i = np.matmul(sym_matrix,G_G[i])
            not_found = True
            j=0
            while not_found:
                if j>=len(G_G):
                    print('Error')
                if g_i[0]==G_G[j,0] and g_i[1]==G_G[j,1] and g_i[2]==G_G[j,2]:
                    RC_G[i] = C_G[k,j]
                    not_found = False
                else: j+=1
                if j==len(G_G) and not_found:
                    print('not found error')
                    RC_G[i] = 0.+0j
                    not_found = False
                    print('g[{}]={}'.format(i,g_i))
        char += np.vdot(RC_G,C_G[k])
    print('charachter = ',char)
    charachters.append([char, sym_name])
