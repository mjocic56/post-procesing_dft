import numpy as np


def char_coef(g_ired,g_red,Nk):
    prod=0
    for i in range(len(Nk)):
        prod+=g_ired[i]*g_red[i]*Nk[i]
    return prod

# Single Group irreps:


char_g1 = np.array([1,1,1,1,1,1,1,1,1,1])
char_g2 = np.array([1,1,-1,-1,1,1,1,-1,-1,1])
char_g12= np.array([2,2,0,0,-1,2,2,0,0,-1])
char_g15= np.array([3,-1,1,-1,0,-3,1,-1,1,0])
char_g25= np.array([3,-1,-1,1,0,-3,1,1,-1,0])

char_g1p= np.array([1, 1, 1, 1, 1,-1,-1,-1,-1,-1])
char_g2p= np.array([1, 1,-1,-1, 1,-1,-1, 1, 1,-1])
char_g12p=np.array([2, 2, 0, 0,-1,-2,-2, 0, 0, 1])
char_g15p=np.array([3,-1, 1,-1, 0, 3,-1, 1,-1, 0])
char_g25p=np.array([3,-1,-1, 1, 0, 3,-1,-1, 1, 0])

Oh_sym = [[char_g1,"$\Gamma_1$"],[char_g2,"$\Gamma_2$"],
          [char_g12,"$\Gamma_{12}$"],
          [char_g15,"$\Gamma_{15}$"],
          [char_g25,"$\Gamma_{25}$"],
          [char_g1p,"$\Gamma^'_1$"],[char_g2p,"$\Gamma^'_2$"],
          [char_g12p,"$\Gamma^'_{12}$"],
          [char_g15p,"$\Gamma^'_{15}$"],
          [char_g25p,"$\Gamma^'_{25}$"]]

Nk = [1,3,6,6,8,1,3,6,6,8]

def Check_irreps(Deg_char,Nk):  
    ir_rep=[]
    for i in range(len(Deg_char)):
        for char_g in Oh_sym:
            c=char_coef(char_g[0],Deg_char[i],Nk)/48
            if (c)==1:
                  ir_rep.append(char_g[1])

    for i in range(len(ir_rep)):
        print("{} irrep: {}".format(ir_rep[i],Deg_char[i]))
    return ir_rep



# Double Group irreps:

char_g6p = np.array([2,0,1,0,1.414,-2,-1,-1.414,2,0,1,0,1.414,-2,-1,-1.414])
char_g7p = np.array([2,0,1,0,-1.414,-2,-1,1.414,2,0,1,0,-1.414,-2,-1,1.414])
char_g6m= np.array([2,0,1,0,1.414,-2,-1,-1.414,-2,0,-1,0,-1.414,2,1,1.414])
char_g7m= np.array([2,0,1,0,-1.414,-2,-1,1.414,-2,0,-1,0,1.414,2,1,-1.414])

char_g8p= np.array([4,0,-1,0,0,-4,1,0,4,0,-1,0,0,-4,1,0])
char_g8m= np.array([4,0,-1,0,0,-4,1,0,-4,0,1,0,0,4,-1,0])

Nk_double = [1,6,8,12,6,1,8,6,1,6,8,12,6,1,8,6]

Oh_sym_double = [[char_g6m,"$\Gamma_6^-$"],[char_g7m,"$\Gamma_7^-$"],
          [char_g8m,"$\Gamma_8^-$"],
          [char_g6p,"$\Gamma_6^+$"],[char_g7p,"$\Gamma_7^+$"],
          [char_g8p,"$\Gamma_8^+$"]]


def Check_irreps_double(Deg_char,Nk_double):  
    ir_rep=[]
    for i in range(len(Deg_char)):
        for char_g in Oh_sym_double:
            c=char_coef(char_g[0],Deg_char[i],Nk_double)/96
            if abs(c-1.0)<1e-3:
                ir_rep.append(char_g[1])
            
    for i in range(len(ir_rep)):
        print("{} irrep: {}".format(ir_rep[i],Deg_char[i]))
    return ir_rep

