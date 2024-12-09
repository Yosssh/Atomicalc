import numpy as np

def TGI(n1,n2):
    numbers = [n1]
    n = n1
    while n < n2:
        n += 1
        numbers.append(n)
    
    return numbers

class Laguerre:
    def __init__(self,p):
        pass #L_p(x)=e^x (d^p)/(dx^p)*(e^(-x)*x^p)

class Asolaguerre:#Tambien se puede usar sum(s=0 to p)[(-1)^s * x^s((p+k)!)^2/((p-s)!*(k+s)!*s!)]
    def __init__(self):#en particular: L_p^k=L_p y L_0^k=k!
        pass #L_p^k(x) = -1^k (d^k)/(dx^k)*L_p+k(x)

class Bra:
    def __init__(self, el):
        super().__init__()
        self.val = np.array(el, dtype=complex)

    def __str__(self):
        return f"<{self.val}|"
    
    def __mul__(self,other):
        if isinstance(other, Ket):
            return np.dot(self.val, other.val)
        elif isinstance(other, Operator):
            return Bra(np.dot(other.M, self.val))
        else:
            print('error producto bra')
    
    def compl_conj(self):#falta añadir que sea el conjugado xD
        return Ket(list(self.val.conj()))

class Ket:
    def __init__(self, el):
        super().__init__()
        self.val = np.zeros(shape=[len(el),1], dtype=complex)
        self.val[:,0] = el

    def __str__(self):
        return f"|{self.val.T[0]}>"

    def compl_conj(self):#falta añadir que sea el conjugado xD
        return Bra(list(self.val.conj().T[0]))
    
    def diad_prod(self, otroket):#toma un conjunto n de numeros y los añade al ket
        pass

class Operator:
    def __init__(self, M):
        self.M = M

    def __mul__(self, other):
        if isinstance(other, Ket):
            return Ket(np.dot(self.M, other.val))

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

class N:
    def __init__(self, n=1):
        self.n = n

    def __str__(self):
        return(self.val)

class L:
    def __init__(self, Qn=1, val = 0):
        self.MAX = Qn-1
        self.rang = TGI(0,self.MAX)
        if val in self.rang:
            self.val = val
        else:
            self.val = 0

    def __str__(self):
        return(self.val)

class S:
    def __init__(self, val=1/2):
        self.val = val

    def __str__(self):
        return(self.val)

class M:
    def __init__(self, MAX=0, val = 0):
        self.MAX = MAX
        self.rang = TGI(-MAX,MAX)
        if val in self.rang:
            self.val = val
        else:
            self.val = 0

    def __str__(self):
        return(self.val)

class OrtBass:
    def __init__(self, basis): #Yo quiero darle algo tipo [a,b,c] y que me instancie los kets
        self.basis = []
        for i in range(len(basis)):
            vec = list(np.zeros(len(basis)))
            vec[i] = 1
            self.basis.append(Ket(vec))

    def __str__(self):
        return str([str(v) for v in self.basis])


class State:
    def __init__(self, basis, coefs):
        self.basis = basis
        self.ket = Ket(coefs)
        self.Proyectors = [Operator((np.outer(ket.val, ket.val.T))) for ket in self.basis]

    def __str__(self):
        letras = ['a', 'b', 'c', 'd', 'e', 'f']
        ch = ''
        for i in range(len(self.basis)):
            if i == 0:
                ch = f'{self.ket.val[i,0]}|a> '
            else:
                ch += f"+ {self.ket.val[i,0]}|{letras[i]}> "
        return ch

#Yo quiero darle a State los números cuanticos y que lo haga todo, cómo?

