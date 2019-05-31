import numpy as np
import re as re

from Symmetry_class import *
from post_procesing_tools import *

def Gather_C_so(path,k_point,details=False,form='stack'):    
    FILE_up = open (path+'/wfc{}_1.xml'.format(k_point),"r")
    FILE_dn = open (path+'/wfc{}_2.xml'.format(k_point),"r")
    nbnd=0
    igwx=0
    WF_U=[]
    WF_D=[]
    if details:
        print("Now performing for k-point={}:\n".format(k_point))
    for line in FILE_up:
        if 'INFO' in line:
            found_nbnd = re.search('nbnd="(.+?)"',line)
            if found_nbnd:
                nbnd = int(found_nbnd.group(1))
        if 'igwx' in line:
            found_igwx=re.search('igwx="(.+?)"',line)
            if found_igwx:
                igwx = int(found_igwx.group(1))
        for n in range(nbnd):
            wf_n=[]
            if '<evc.{} type='.format(n+1) in line:
                for i in range(igwx):
                    nextLine=next(FILE_up)
                    A = (nextLine.replace("\n","")).split(",")
                    C = complex(float(A[0]),float(A[1]))

                    wf_n.append(C)
                WF_U.append(wf_n)
    for line in FILE_dn:
        if 'Info' in line:
            found_nbnd = re.search('nbnd="(.+?)"',line)
            if found_nbnd:
                nbnd = int(found_nbnd.group(1))
        if 'igwx' in line:
            found_igwx=re.search('igwx="(.+?)"',line)
            if found_igwx:
                igwx = int(found_igwx.group(1))
        for n in range(nbnd):
            wf_n=[]
            if '<evc.{} type='.format(n+1) in line:
                for i in range(igwx):
                    nextLine=next(FILE_dn)
                    A = (nextLine.replace("\n","")).split(",")
                    C = complex(float(A[0]),float(A[1]))

                    wf_n.append(C)
                WF_D.append(wf_n)
    if details:
        print("\tReading wave-fuctions: DONE\n")
    WF_U=np.array(WF_U)
    WF_D=np.array(WF_D)
    FILE_up.close()
    FILE_dn.close()
    if form == 'stack':
        return np.array([WF_U,WF_D])
    if form == 'interweave':
        shape0 = WF_U.shape[0] + WF_D.shape[0]
        shape1 = WF_U.shape[1]
        WF = np.empty((shape0,shape1), dtype=WF_U.dtype)
        WF[0::2] = WF_U; WF[1::2] = WF_D
        return WF
    if form == 'rack':
        shape0 = WF_U.shape[0]
        shape1 = WF_U.shape[1]
        WF = np.empty((2*shape0,shape1), dtype=WF_U.dtype)
        WF[0:shape0] = WF_U; WF[shape0:2*shape0] = WF_D
        return WF 

def Gather_G_so(path,k_point,details=False):
    FILE = open (path+'/grid{}.xml'.format(k_point),"r")
    size=0
    GV=[]
    if details:
        print("Gathering G-vectors for k-point={}:\n".format(k_point))
    for line in FILE:
        if '<GRID type=' in line:
            size = re.search('size="(.+?)"',line)
            if size:
                size = int(int(size.group(1))/3)
            for i in range(size):
                nextLine=next(FILE)
                g = (nextLine.replace("\n","")).split()
                G_i = np.array([int(g[0]),int(g[1]),int(g[2])])
                GV.append(G_i)
    if details:
        print("\tReading G-vectors: DONE\n")
    GV=np.array(GV)
    return GV

def Gather_Degeneracies_so(E_list,functions=False):
    deg_e=[]
    deg_u=[]
    index=0
    while index+4<E_list.shape[0]:
        if abs(E_list[index]-E_list[index+1])<1e-4:
            index += 1
            if abs(E_list[index]-E_list[index+1])<1e-4:
                index += 1
                if abs(E_list[index]-E_list[index+1])<1e-4:
                    index += 1
                    deg_e.append([E_list[index],4,'index: {}-{}'.format(index-3,index)])
                    deg_u.append([index-3,index-2,index-1,index])
                else:
                    deg_e.append([E_list[index],3,'index: {}-{}'.format(index-2,index)])
                    deg_u.append([index-2,index-1,index])
            else:
                deg_e.append([E_list[index],2,'index: {}-{}'.format(index-1,index)])
                deg_u.append([index-1,index])
        else:
            deg_e.append([E_list[index],1,'index: {}'.format(index)])
            deg_u.append([index])
        index +=1
    if index == E_list.shape[0] - 2:
        deg_e.append([E_list[index],2,'index: {}-{}'.format(index,index+1)])
        deg_u.append([index,index+1])
    if functions==False:
        return np.array(deg_e)
    else:
        return (deg_u),np.array(deg_e)
    
def Gather_Symmetries_double(path,file_name="index.xml"):
    full_path=path+"/"+file_name
    FILE = open (full_path,"r")
    n_sym=0
    inv_sym = False
    sym_list=[]
    sym_list_d=[]
    sigma_list=[]
    sigma_list_d=[]
    for line in FILE:
        if '<symmops' in line:
            n_sym = re.search('nsym="(.+?)"',line)
            if n_sym:
                n_sym = int(n_sym.group(1))
    for i in range(1,n_sym+1):
        FILE = open (full_path,"r")
        for line in FILE:
            if '<info.{} name'.format(i) in line:
                #Single-group
                sym_name = re.search('name="(.+?)"',line)
                sym_name = str(sym_name.group(1))
                nextLine=next(FILE)
                mat_type = re.search('type="integer"',nextLine)
                nextLine=next(FILE)
                r = nextLine.split()
                R1 = [int(r[0]),int(r[1]),int(r[2])]
                nextLine=next(FILE)
                r = nextLine.split()
                R2 = [int(r[0]),int(r[1]),int(r[2])]
                nextLine=next(FILE)
                r = nextLine.split()
                R3 = [int(r[0]),int(r[1]),int(r[2])]
                mat = np.array([R1,R2,R3])
                R = Symmetry(mat,sym_name)
                Rd= Symmetry(mat,sym_name+"+d")
                sym_list.append(R)
                sym_list_d.append(Rd)
                #Double-group:
                nextLine=next(FILE)
                mat_type = re.search('type="complex"',nextLine)
                nextLine=next(FILE)
                s = nextLine.split()
                S1 = [eval(s[0],{'sqrt':np.sqrt}),eval(s[1],{'sqrt':np.sqrt})]
                nextLine=next(FILE)
                s = nextLine.split()
                S2 = [eval(s[0],{'sqrt':np.sqrt}),eval(s[1],{'sqrt':np.sqrt})]
                sigma = np.array([S1,S2],dtype=complex)
                S = Symmetry(sigma,sym_name)
                Sd = Symmetry(-sigma,sym_name+"+d")
                sigma_list.append(S)
                sigma_list_d.append(Sd)
    sym_list.extend(sym_list_d)
    sigma_list.extend(sigma_list_d)
    return sym_list, sigma_list

def Gather_Charachters_double(k_start,k_stop,symmetry_list,sigma_list,G_R,C_Rup,C_Rdn,k_Rpoint,details=False):
    charachters=[]
    dim_G = 2*abs(max(np.amax(G_R),np.amin(G_R),key=abs))+1
    #dim_G = int(dim_G)
    if details: print('Dimension of 3x3 square matrix is:', dim_G)
    M=np.zeros((dim_G,dim_G,dim_G),dtype=int)
    if details: print('Creating G-matrix')
    for i in range(len(G_R)):
        n1,n2,n3 = G_R[i,0],G_R[i,1],G_R[i,2]
        M[n1,n2,n3] = i
    if details: print('G-matrix: Done')
    for s in [0,1,5,12,14,48,52,62,24,25,28,36,38,72,76,87]:
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
        char=0
        for k in range(k_start,k_stop+1):
            RC_Rup = np.zeros(len(G_R),dtype=complex)
            RC_Rdn = np.zeros(len(G_R),dtype=complex)
            
            for i in range(len(G_R)):
                g_i = np.matmul(sym_matrix,G_R[i]+k_Rpoint)-k_Rpoint
                n1,n2,n3=int(g_i[0]),int(g_i[1]),int(g_i[2])
                index = M[n1,n2,n3]
                RC_Rup[i] = C_Rup[k,index]
                RC_Rdn[i] = C_Rdn[k,index]
            
            up_up = np.vdot(RC_Rup,C_Rup[k])*aa
            up_dn = np.vdot(RC_Rup,C_Rdn[k])*ab
            dn_up = np.vdot(RC_Rdn,C_Rup[k])*ba
            dn_dn = np.vdot(RC_Rdn,C_Rdn[k])*bb
            
            char += up_up + up_dn + dn_up + dn_dn
        if details: print('charachter = ',round(char,3))
        charachters.append([round(char.real,3), sym_name])
    if details: print('\nDONE.\n')
    return np.array(charachters)

def Check_Orthonormality_spin(WF_up,WF_dn, tol=1e-13,details=False):
    # check orthonormality
    nbnd = len(WF_up)
    Ort_MAT=np.zeros((nbnd,nbnd),dtype=complex)
    for i in range(nbnd):
        for j in range(i,nbnd):
            A_up=np.array(WF_up[i])
            A_dn=np.array(WF_dn[i])
            B_up=np.array(WF_up[j])
            B_dn=np.array(WF_dn[j])
            Ort_MAT[i][j]=np.vdot(A_up,B_up)+np.vdot(A_dn,B_dn)
            Ort_MAT[j][i]=Ort_MAT[i][j]
    if details:
        print("\tVector product is: DONE")
    
    orthonormal=True
    #print("\nDiagonal elements:\n")
    tol_global=tol
    index = [0,0]
    for i in range(nbnd):
        if abs(abs(Ort_MAT[i][i])-1.0)<tol:
            pass#print (1)
        else:
            if details:
                print ("\n\tMatrix element [{},{}] is:{}".format(i,i,Ort_MAT[i][i]))
            orthonormal=False
            ort_new=False
            tol_new=tol
            while ort_new==False:
                    tol_new=tol_new*10
                    if tol_global<tol_new: tol_global=tol_new; index=[i,i]
                    if abs(abs(Ort_MAT[i][i])-1.0)<tol_new:
                        if details:
                            print("\tTolerance for {}".format(tol_new))
                        ort_new=True
    #print("\nOff-diagonal elements:\n")
    for i in range(nbnd):
        for j in range(i+1,nbnd):
            if abs(Ort_MAT[i][j])<tol:
                pass#print(0)
            else:
                if details:
                    print ("\n\tMatrix element: {},{} is:{}".format(i,j,abs(Ort_MAT[i][j])))
                orthonormal=False
                ort_new=False
                tol_new=tol
                while ort_new==False:
                    tol_new=tol_new*10
                    if tol_global<tol_new: tol_global=tol_new; index=[i,j]
                    if abs(Ort_MAT[i][j])<tol_new:
                        if details: 
                            print("\tTolerance for {}".format(tol_new))
                        ort_new=True                    
    print('\nOrthonormality {} for tolerance {}'.format(orthonormal,tol))
    if not orthonormal:
        print('Orthonormal for tolerance {} in index {}\n'.format(tol_global,index)) 