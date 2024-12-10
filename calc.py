import numpy as np
import itertools

def TGI(n1,n2):
    numbers = [n1]
    n = n1
    while n < n2:
        n += 1
        numbers.append(n)
    
    return numbers

def Qprod(uno, otro):
    if uno.val == otro.val:
        return 1
    else:
        return 0

class Laguerre:
    def __init__(self,p):
        pass #L_p(x)=e^x (d^p)/(dx^p)*(e^(-x)*x^p)

class Asolaguerre:#Tambien se puede usar sum(s=0 to p)[(-1)^s * x^s((p+k)!)^2/((p-s)!*(k+s)!*s!)]
    def __init__(self):#en particular: L_p^k=L_p y L_0^k=k!
        pass #L_p^k(x) = -1^k (d^k)/(dx^k)*L_p+k(x)

class HVec:#Aqui estoy definiendo bien los kets con su producto pertinente!!!!!!
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
        self.etiqueta = etiqueta
        self.vec = vec

    def __add__(self,other):
        if isinstance(other, Ket):
            return(self.vec+other.vec)
        else:
            return 'error in ket suma'

    def __str__(self):
        ch =f"|"
        for i,v in enumerate(self.etiqueta):
            if i == 0:
                ch += f"{v}"
            else:
                ch += f",{v}"
        return f"{ch}>"
    
    def to_bra(self):
        return Bra(self.etiqueta, self.vec)

    def diad_prod(self, otroket):#toma un conjunto n de numeros y los añade al ket
        pass

class Operator:
    def __init__(self, M):
        self.M = M

    def __mul__(self, other):
        if isinstance(other, Ket):
            return Ket(np.dot(self.M, other.val))#Hay que cambiar esto porque los kets y bras ahora son mas fumada (ver HVec)

        else:
            print('Error en producto operador')

class StaionarytNormalizedWFunction:
    def __init__(self):#Lo uyo sería separarlo en parte radial y angular
        pass#a_0^(-3/2)*(2)/(n**2)sqrt[(n-l-1)!/((n+l)!^3)]((2r)/(na_o))^l * exp(-r/na_0) * L_n-l-1^2l+1((2r)/(na_0))Y_l^m

class Orbital:
    def __init__(self, n, l, m, rad_wf=None):
        self.n = n
        self.l = l
        self.m = m

class Atom:
    def __init__(self, Z, N):
        self.q = Z
        self.n = N

class OrtBass:
    def __init__(self, n=[1], l=[], ml=[], s=[], ms=[]): #Yo quiero darle los números cuánticos y tener todos los posibles estados sin tener en cuenta restricciones
        self.base = [] #lista de kets de la base
        self.n = n
        if len(l) == 0:
            self.l = TGI(0,n[0]-1)
        else:
            self.l = l
        if len(ml) == 0:
            self.ml = TGI(-max(self.l),max(self.l))
        else:
            self.ml = ml
        auxiliar = list(itertools.product(*[self.n, self.l, self.ml]))
        base = [comb for comb in auxiliar if comb[1]>=abs(comb[2]) and comb[0]==self.n[0]] #etiquetas de números
        base_vecs = [np.zeros(len(base)) for _ in base] #arrays tipo [1,0,0]
        for i,b in enumerate(base_vecs):
            b[i] = 1

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
        pass



#Yo quiero darle a State los números cuanticos y que lo haga todo, cómo?

base = OrtBass([2]) #Creamos una base de vectores para lo números cuánticos dados (los especificados se mantienen fijos, los que no se tomarán todos los valores válidos)
print(base)
estado = State(base)
estado.nicks = True
print(estado)
