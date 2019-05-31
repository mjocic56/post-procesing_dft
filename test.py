import numpy as np
import re as re

from Symmetry_class import *
from post_procesing_tools import *

def Gather_irrep2(path,file_name="index.xml"):
    full_path=path+"/"+file_name
    FILE = open (full_path,"r")
    n_sym=0
    inv_sym = False
    sigma_list=[]
    for line in FILE:
        if '<symmops' in line:
            n_sym = re.search('nsym="(.+?)"',line)
            if n_sym:
                n_sym = int(n_sym.group(1))
    for i in range(1,n_sym+1):
        FILE = open (full_path,"r")
        for line in FILE:
            if '<info.{} name'.format(i) in line:
                #Double-group:
                sym_name = re.search('name="(.+?)"',line)
                sym_name = str(sym_name.group(1))
                #print(sym_name)
                nextLine=next(FILE)
                mat_type = re.search('type="complex"',nextLine)
                nextLine=next(FILE)
                s = nextLine.split()
                #print(s)
                S1 = [eval(s[0],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[1],{'exp':np.exp,'pi':np.pi})]
                nextLine=next(FILE)
                s = nextLine.split()
                S2 = [eval(s[0],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[1],{'exp':np.exp,'pi':np.pi})]
                sigma = np.array([S1,S2],dtype=complex)
                S = Symmetry(sigma,sym_name)
                sigma_list.append(S)
    return sigma_list
def Gather_irrep4(path,file_name="index.xml"):
    full_path=path+"/"+file_name
    FILE = open (full_path,"r")
    n_sym=0
    inv_sym = False
    sigma_list=[]
    for line in FILE:
        if '<symmops' in line:
            n_sym = re.search('nsym="(.+?)"',line)
            if n_sym:
                n_sym = int(n_sym.group(1))
    for i in range(1,n_sym+1):
        FILE = open (full_path,"r")
        for line in FILE:
            if '<info.{} name'.format(i) in line:
                #Double-group:
                sym_name = re.search('name="(.+?)"',line)
                sym_name = str(sym_name.group(1))
                #print(sym_name)
                nextLine=next(FILE)
                mat_type = re.search('type="complex"',nextLine)
                nextLine=next(FILE)
                s = nextLine.split()
                #print(s)
                S1 = [eval(s[0],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[1],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[2],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[3],{'exp':np.exp,'pi':np.pi})]
                #print(S1)
                nextLine=next(FILE)
                s = nextLine.split()
                S2 = [eval(s[0],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[1],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[2],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[3],{'exp':np.exp,'pi':np.pi})]
                nextLine=next(FILE)
                s = nextLine.split()
                S3 = [eval(s[0],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[1],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[2],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[3],{'exp':np.exp,'pi':np.pi})]
                nextLine=next(FILE)
                s = nextLine.split()
                S4 = [eval(s[0],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[1],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[2],{'exp':np.exp,'pi':np.pi}),\
                      eval(s[3],{'exp':np.exp,'pi':np.pi})]
                sigma = np.array([S1,S2,S3,S4],dtype=complex)
                S = Symmetry(sigma,sym_name)
                sigma_list.append(S)
    return sigma_list
    

def Gather_Matrices_double(k_start,k_stop,
                           symmetry_list,sigma_list,
                           G_R,C_Rup,C_Rdn,k_Rpoint,
                           details=False,all_elements=False,
                          elements_list=0):
    Matrices=[]
    dim_G = 2*abs(max(np.amax(G_R),np.amin(G_R),key=abs))+1
    #dim_G = int(dim_G)
    if details: print('Dimension of 3x3x3 square matrix is:', dim_G)
    M=np.zeros((dim_G,dim_G,dim_G),dtype=int)
    if details: print('Creating G-matrix')
    for i in range(len(G_R)):
        n1,n2,n3 = G_R[i,0],G_R[i,1],G_R[i,2]
        M[n1,n2,n3] = i
    if details: print('G-matrix: Done')
    if all_elements:
        s_list = range(len(symmetry_list))
    else:
        s_list = [0,1,4,12,14,48,52,62,24,25,28,36,38,72,76,87]
    if elements_list!=0:
        s_list = elements_list
    for s in s_list:
        sym_name=symmetry_list[s].operation
        sym_matrix=symmetry_list[s].matrix
        sigma_matrix=sigma_list[s].matrix
        if details: 
            print('\nNow performing for '+ sym_name)
            print('\nSymm.Op. matrix:\n{}\nSpin.Op. matrix:\n{}'.format(sym_matrix,sigma_matrix))
        aa = sigma_matrix[0,0]
        ab = sigma_matrix[0,1]
        ba = sigma_matrix[1,0]
        bb = sigma_matrix[1,1]
        
        dim = k_stop - k_start +1
        full_mat = np.zeros((dim,dim),dtype=complex)
        for k in range(k_start,k_stop+1):
            RC_Rup = np.zeros(len(G_R),dtype=complex)
            RC_Rdn = np.zeros(len(G_R),dtype=complex)
            
            for i in range(len(G_R)):
                g_i = np.matmul(sym_matrix,G_R[i]+k_Rpoint)-k_Rpoint
                n1,n2,n3=int(g_i[0]),int(g_i[1]),int(g_i[2])
                index = M[n1,n2,n3]
                RC_Rup[i] = C_Rup[k,index]
                RC_Rdn[i] = C_Rdn[k,index]
                for n in range(k_start,k_stop+1):
                    
                    up_up = np.vdot(RC_Rup,C_Rup[n])*aa
                    up_dn = np.vdot(RC_Rup,C_Rdn[n])*ab
                    dn_up = np.vdot(RC_Rdn,C_Rup[n])*ba
                    dn_dn = np.vdot(RC_Rdn,C_Rdn[n])*bb
                
                    mat_comp = up_up + up_dn + dn_up + dn_dn
                    if abs(mat_comp.real)<1e-9:
                        mat_comp = mat_comp - mat_comp.real
                    if abs(mat_comp.imag)<1e-9:
                        mat_comp = mat_comp -1j*mat_comp.imag
                    
                    full_mat[k-k_start,n-k_start] = mat_comp
        if details: print('Matrix = \n',full_mat)
        Matrices.append(Symmetry(full_mat,sym_name))
    if details: print('\nDONE.\n')
    return np.array(Matrices)


def Compute_Perturbation_spin(G_vec,C_terms_up,C_terms_dn,
                              ext_basis_up,ext_basis_dn,
                              E_terms,E_ext,kkr,m,n):
        new_basis_up = np.vstack((ext_basis_up,C_terms_up))
        new_basis_dn = np.vstack((ext_basis_dn,C_terms_dn))
        dim = len(ext_basis_up)
        cgc = 0
        Eavg = (E_terms[m]+E_terms[n])/2
        for l in range(len(E_ext)):
            delta_E1 =1/(Eavg - E_ext[l])
            ml_up = np.dot(kkr,Compute_Matrix_Element(G_vec,new_basis_up,dim+m,l))
            ml_dn = np.dot(kkr,Compute_Matrix_Element(G_vec,new_basis_dn,dim+m,l))
            ln_up = np.dot(kkr,Compute_Matrix_Element(G_vec,new_basis_up,l,dim+n))
            ln_dn = np.dot(kkr,Compute_Matrix_Element(G_vec,new_basis_dn,l,dim+n))
            kpi = ml_up + ml_dn
            kpj = ln_up + ln_dn
            cgc += kpi*kpj*delta_E1
        return cgc
    
def Compute_KdP2_spin(k_mesh,k_point,E_k,E_ext,G_space,
                      C_terms_up,C_terms_dn,
                      ext_basis_up,ext_basis_dn,
                      a,units='Hartree',details=False):
    uni=['Hartree',1.0]
    if units=='Rydberg':
        uni=['Rydberg',2.0]
    if units=='Electronvolt':
        uni=['Electronvolt',27.211396132]
    tpiba = 2*np.pi/a
    tpiba2=tpiba**2
    tpiba4=tpiba2**2
    H=[]
    H_eig=[]
    dim = len(C_terms_up)
    for k_i in range(len(k_mesh)):
        k2 = np.dot(k_mesh[k_i],k_mesh[k_i])
        kkR = k_mesh[k_i]-k_point
        H_i = np.zeros((dim,dim),dtype=complex)
        for m in range(dim):
            for n in range(dim):
                KdP2 = Compute_Perturbation_spin(G_space,C_terms_up,C_terms_dn,
                                                 ext_basis_up,ext_basis_dn,
                                                 E_k,E_ext,kkR,m,n)
                H_i[m][n] = tpiba4*KdP2
        H_i_eig = (np.linalg.eigvalsh(H_i))
        H.append(H_i)
        H_eig.append(H_i_eig)
        if details:
            print("eigenvalues of KdP2_{} in units of {}".format(k_i+1,uni[0]))
            print((H_i_eig)*uni[1])
    print('Calculating KdP2: DONE\n')
    return np.array(H)*uni[1], np.array(H_eig)*uni[1]

def Compute_KdP1_spin(k_mesh,k_point,G_space,
                      C_terms_up,C_terms_dn,a,units='Hartree',details=False):
    uni=['Hartree',1.0]
    if units=='Rydberg':
        uni=['Rydberg',2.0]
    if units=='Electronvolt':
        uni=['Electronvolt',27.211396132]
    tpiba = 2*np.pi/a
    tpiba2=tpiba**2
    tpiba4=tpiba2**2
    H=[]
    H_eig=[]
    dim = len(C_terms_up)
    for k_i in range(len(k_mesh)):
        kkR = k_mesh[k_i]-k_point
        H_i = np.zeros((dim,dim),dtype=complex)
        for m in range(dim):
            for n in range(dim):
                KdP_up = Compute_Matrix_Element(G_space,C_terms_up,m,n)
                KdP_dn = Compute_Matrix_Element(G_space,C_terms_dn,m,n)
                H_i[m][n] = tpiba2*kkR.dot(KdP_up+KdP_dn)
        H_i_eig = (np.linalg.eigvalsh(H_i))
        H.append(H_i)
        H_eig.append(H_i_eig)
        if details:
            print("eigenvalues of KdP_{} in units of {}".format(k_i+1,uni[0]))
            print((H_i_eig)*uni[1])
    print('Calculating KdP1: DONE\n')
    return np.array(H)*uni[1], np.array(H_eig)*uni[1]