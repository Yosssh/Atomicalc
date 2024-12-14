import numpy as np
import itertools
import ajustes

def TGI(n1,n2):
    numbers = [n1]
    n = n1
    while n < n2:
        n += 1
        numbers.append(n)
    
    return numbers

class StaionarytNormalizedWFunction:
    def __init__(self):#Lo uyo sería separarlo en parte radial y angular
        pass#a_0^(-3/2)*(2)/(n**2)sqrt[(n-l-1)!/((n+l)!^3)]((2r)/(na_o))^l * exp(-r/na_0) * L_n-l-1^2l+1((2r)/(na_0))Y_l^m

class Laguerre:
    def __init__(self,p):
        pass #L_p(x)=e^x (d^p)/(dx^p)*(e^(-x)*x^p)

class Asolaguerre:#Tambien se puede usar sum(s=0 to p)[(-1)^s * x^s((p+k)!)^2/((p-s)!*(k+s)!*s!)]
    def __init__(self):#en particular: L_p^k=L_p y L_0^k=k!
        pass #L_p^k(x) = -1^k (d^k)/(dx^k)*L_p+k(x)

class HVec:#Aqui iba a definir bien los kets con su producto pertinente!!!!!!
    def __init__(self, v):
        self.v = v
        self.n = v[0]
        self.n = v[1]
        self.m = v[2]
        if len(v) == 5:
            self.s = v[3]
            self.ms = v[4]
        else:
            pass

        def __mul__(self, other):
            if isinstance(other, HVec):
                if self.v == other.v:
                    return 1
                else:
                    return 0
            
            else:
                return 'Error en HVec __mul__'

class Bra:
    def __init__(self, etiqueta, vec):
        self.etiqueta = etiqueta
        self.vec = vec

    def __add__(self, other):
        if isinstance(other, Bra):
            return (self.vec+other.vec)
        else:
            return 'error en bra suma'

    def __str__(self):
        ch =f"<"
        for i,v in enumerate(self.etiqueta):
            if i == 0:
                ch += f"{v}"
            else:
                ch += f",{v}"
        return f"{ch}|"
    
    def __mul__(self,other):
        if isinstance(other, Ket):
            return np.dot(self.vec, other.vec)
        elif isinstance(other, Operator):
            return Bra(np.dot(other.M, self.val))
        else:
            print('error producto bra')
    
    def to_ket(self):
        return(Ket(self.etiqueta, self.vec))

class Ket:
    def __init__(self, etiqueta, vec):
        self.etiqueta = np.array(etiqueta)
        self.vec = vec
        self.showing = [True for _ in self.etiqueta]

    def __add__(self,other):
        if isinstance(other, Ket):
            return(self.vec+other.vec)
        else:
            return 'error in ket suma'

    def __str__(self):
        ch =f"|"
        visibletiquets = [e for i,e in enumerate(self.etiqueta) if self.showing[i]]            
        for i,v in enumerate(visibletiquets):
            if i == 0:
                ch += f"{v}"
            else:
                ch += f",{v}"
        return f"{ch}>"
    
    def to_bra(self):
        return Bra(self.etiqueta, self.vec)

    def __mull__(self,other):
        if isinstance(other, float) or isinstance(other,int):
            return self.vec*other

    def actualizar_ket(self, etiqueta, vec):
        return Ket(etiqueta,vec)

    def diad_prod(self, other):#Devuelve unos kets cuyo .vec no esta definido
        if isinstance(other, Ket):
            new_etiqueta = self.etiqueta+other.etiqueta
            return Ket(new_etiqueta,None)

class Orbital:
    def __init__(self, n, l, occ_num):
        self.n = n
        self.l = l
        self.occ_num = occ_num

    def __str__(self):
        return (f"Orbital {self.n}{ajustes.nomenclatura_orbitalesinv.get(self.l)}")

    def is_closed(self):
        if self.occ_num == 2*(2*self.l + 1):
            return True 
        else:
            return False

class Atom:
    def __init__(self, elec_conf):
        self.numbers = self.get_elec_conf(elec_conf)
        self.orbitals = [Orbital(n[0],n[1],n[2]) for n in self.numbers]
        self.open_orbs = [o for o in self.orbitals if not o.is_closed()]

    def __str__(self):
        ch = f"------Orbitales-----------------"
        for i,o in enumerate(self.orbitals):
            ch += f'\n{str(o)}      nº ocupación: {self.numbers[i][2]}'
        ch += "\n-------------------------------"
        return ch
        

    def get_elec_conf(self, string): #Pedimos que la configuracion electronica se introduzca separada por espacios ej: 1s2 2p
        orbs = string.split(' ')
        orbs_numbers = []
        for o in orbs:
            if len(o) == 2:
                orbs_numbers.append([int(o[0]),ajustes.nomenclatura_orbitales.get(o[1]),1])
            else:
                orbs_numbers.append([int(o[0]), ajustes.nomenclatura_orbitales.get(o[1]), int(o[2])])

        return(orbs_numbers)
    
    def get_open_nums(self):
        n = list(set([o.n for o in self.open_orbs]))
        l = list(set([o.l for o in self.open_orbs]))
        return(n,l)
    
    def get_base(self,spin=False):#Para estudiar los electrones de las capas abiertas
        n,l = self.get_open_nums()
        if not spin:
            return OrtBass(n,l)
        if spin:
            return OrtBass(n,l,s=[0.5])

class OrtBass:
    def __init__(self, n=[1], l=[], ml=[], s=[], ms=[]): #Yo quiero darle los números cuánticos y tener todos los posibles estados sin tener en cuenta restricciones
        self.base = [] #lista de kets de la base
        self.n = n
        self.spin = None
        if len(l) == 0:
            self.l = TGI(0,n[0]-1)
        else:
            self.l = l
        if len(ml) == 0:
            self.ml = TGI(-max(self.l),max(self.l))
        else:
            self.ml = ml
        if len(s) == 0 and len(ms) == 0:
            self.spin = False
            auxiliar = list(itertools.product(*[self.n, self.l, self.ml]))
            base = [comb for comb in auxiliar if comb[1]>=abs(comb[2]) and comb[1]<comb[0]] #etiquetas de números
            base_vecs = [np.zeros(len(base)) for _ in base] #arrays tipo [0,0,0]
            for i,b in enumerate(base_vecs):
                b[i] = 1 #Se hace 1 un único elemento
        else:
            self.spin = True
            self.s = s
            if len(ms) == 0:
                self.ms = TGI(-s[0],s[0])
            else:
                self.ms = ms
            auxiliar = list(itertools.product(*[self.n, self.l, self.ml,self.s,self.ms]))
            base = [comb for comb in auxiliar if comb[1]>=abs(comb[2]) and comb[1]<comb[0]] #etiquetas de números
            base_vecs = [np.zeros(len(base)) for _ in base] #arrays tipo [0,0,0]
            for i,b in enumerate(base_vecs):
                b[i] = 1 #Se hace 1 un único elemento

        for etiqueta,vec in zip(base,base_vecs):
            self.base.append(Ket(etiqueta,vec))

    def __str__(self):
        ch = ''
        for k in self.base:
            ch += f"{k}\n"
        return ch

class State: #Creo que hay que darle una vuelta a esta o la clase Ket para unificarlas y reducir espacio (además de que por definicion un estado es un ket)
    def __init__(self, basis, coefs=None):
        self.basis = basis
        if coefs == None:
            self.coefs = np.ones(len(basis.base))
        else: 
            self.coefs = coefs
        self.nicks = False
        self.Proyectors = [Operator((np.outer(ket.vec, ket.to_bra().vec))) for ket in self.basis.base]

    def __str__(self):
        if self.nicks:
            letras = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u']
            ch = ''
            for i in range(len(self.basis.base)):
                if i == 0:
                    ch = f'{self.coefs[i]}|a> '
                else:
                    ch += f"+ {self.coefs[i]}|{letras[i]}> "
            return ch
        else:
            ch = ''
            for i,k in enumerate(self.basis.base):
                if i == 0:
                    ch += f'{self.coefs[i]}{k} '
                else:
                    ch += f"+ {self.coefs[i]}{k} "
            return ch

    def Slater_Det(self):
        combs = []
        for k1 in self.basis.base:
            for k2 in self.basis.base:
                combs.append(k1.diad_prod(k2))

        combs_posibles = list(set([c for c in combs if not Pauli_exclusion(c)]))

        return combs_posibles

class Operator:
    def __init__(self, transformation, base, M=None, fun=None):
        self.t = transformation
        self.fun = fun
        if M != None:
            self.M = M
        else:
            self.M = self.get_M(base)


    def __mul__(self, other):
        if isinstance(other, Ket):
            return Ket(other.etiqueta+self.t,np.dot(self.M, other.vec))

        else:
            print('Error en producto operador')
    
    def get_M(self, base):
        kets_base = base.base
        kets_base_e = [list(k.etiqueta) for k in kets_base]
        M = []
        for ki in kets_base:
            fila = []
            for kj in kets_base:
                mij = self.fun(kj)
                etiqueta = kj.etiqueta + self.t
                if list(etiqueta) in kets_base_e and list(etiqueta)==list(ki.etiqueta):
                    fila.append(mij)
                else:
                    fila.append(0)
            M.append(fila)
        return np.array(M)

def Pauli_exclusion(ket):#esto es un ket de 10 numeros (los dos electrones)
    if ket.etiqueta[0]==ket.etiqueta[5] and ket.etiqueta[1]==ket.etiqueta[6] and ket.etiqueta[2]==ket.etiqueta[7] and ket.etiqueta[4]==ket.etiqueta[9]:
        return True
    else:
        return False



##########################################################################################
caso = Atom('2p')

base = caso.get_base(spin=True)
base_kets = base.base
##########################################################################################
def Lsqr_exp(ket):
    c = ket.etiqueta[1]*(ket.etiqueta[1]+1)
    return c
lcuad = Operator(np.array([0,0,0,0,0]), base, fun=Lsqr_exp)

def Lplus_exp(ket):
    c = np.sqrt(ket.etiqueta[1]*(ket.etiqueta[1]+1) - ket.etiqueta[2]*(ket.etiqueta[2]+1))
    return c
lplus = Operator(np.array([0,0,1,0,0]), base, fun=Lplus_exp)


def Lminus_exp(ket):
    c = np.sqrt(ket.etiqueta[1]*(ket.etiqueta[1]+1) - ket.etiqueta[2]*(ket.etiqueta[2]-1))
    return c
lminus = Operator(np.array([0,0,-1,0,0]), base, fun=Lminus_exp)


def Lz_exp(ket):
    c = ket.etiqueta[2]
    return c
lz = Operator(np.array([0,0,0,0,0]), base, fun=Lz_exp)

LxM = lplus.M/2 + lminus.M/2  
LyM = 1j*lminus.M/2 - 1j*lplus.M/2
#otroL2 = np.dot(lz.M,lz.M) + np.dot(LxM,LxM) + np.dot(LyM,LyM)#SALEEEEEE!!!!
##########################################################################################

def Ssqr_exp(ket):
    try:
        c = ket.etiqueta[3]*(ket.etiqueta[3]+1)
        return c
    except:
        print("XDDDD")
Scuad = Operator([0,0,0,0,0], base, fun=Ssqr_exp)
def Splus_exp(ket):
    try:
        c = np.sqrt(ket.etiqueta[3]*(ket.etiqueta[3]+1) - ket.etiqueta[4]*(ket.etiqueta[4]+1))
        return c
    except:
        print("XDDDD")
Splus = Operator([0,0,0,0,1], base, fun=Splus_exp)

def Sminus_exp(ket):
    try:
        c = np.sqrt(ket.etiqueta[3]*(ket.etiqueta[3]+1) - ket.etiqueta[4]*(ket.etiqueta[4]-1))
        return c
    except:
        print("XDDDD")
Sminus = Operator([0,0,0,0,-1], base, fun=Sminus_exp)

def Sz_exp(ket):
    try:
        c = ket.etiqueta[4]
        return c
    except:
        print("XDDDD")
Sz = Operator([0,0,0,0,0], base, fun=Sz_exp)
#######################################
print(f"base:\n{base}")
print(f"-----------\nL2:\n{lcuad.M}")
print(f"-----------\nL-:\n{lminus.M}")
print(f"-----------\nL+:\n{lplus.M}")
print(f"-----------\nLz:\n{lz.M}")
print(f"-----------\nS2:\n{Scuad.M}")
print(f"-----------\nS-:\n{Sminus.M}")
print(f"-----------\nS+:\n{Splus.M}")
print(f"-----------\nSz:\n{Sz.M}")
#######################################



miket = base_kets[4]
miket.showing = [False,True,True,False,True]
print(f"\nVamos a ver cómo actuan todos estos operadores sobre {miket}")

print()