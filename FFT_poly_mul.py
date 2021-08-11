import numpy as np

def pretty_str(complex_list: list):
    s="("+str(list(np.round(complex_list, 4)))+')'
    s=s.replace('(','').replace(')','')
    s = s.replace('[', '( ').replace(']', ')')
    s=s.replace(' +0',' 0').replace(' -0',' 0')
    s = s.replace(' -0j', ' 0').replace(' +0j', ' 0').replace(' 0j', ' 0')
    s = s.replace('-0j', '').replace('+0j', '')
    s = s.replace('j', 'i')
    return s

def print_polynom(c: list, polinom_name: str):
    n=len(c)
    print(polinom_name+'(x)=',end='')
    for i in range(n-1,-1,-1):
        ci_s = ("(" + str(np.round(c[i], 4)) + ')').replace('(', '').replace(')', '').replace(' +0', ' 0').replace(' -0', ' 0').replace(' -0j', ' 0').replace(' +0j', ' 0').replace(' 0j', ' 0').replace('-0j', '').replace('+0j', '').replace('j', 'i')
        print(ci_s+('x^'+str(i) if i>0 else ''),end=('\n' if i==0 else '+'))


def abFixLen_and_get_n(a: list, b: list):
    # n=deg of c round up to 2**q (q is some natural)
    deg = lambda v: len(v) - 1
    deg_c = deg(a) + deg(b)
    n = 1
    while n < deg_c:
        n *= 2

    # padding a and b from right
    for i in range(n - len(a)):
        a.append(0)
    for i in range(n - len(b)):
        b.append(0)

    return n

def emptyList(size: int):
    return [0 for i in range(size)]

#Recursive FFT
def FFT(n: int, a: list):
    if n==1:
        return [a[0]]

    t = FFT(n//2,[a[i] for i in range(0,n,2)])
    u = FFT(n//2,[a[i] for i in range(1,n,2)])

    e = np.e; pi = np.pi
    W = 1
    W_n = e**( 1j*((2*pi)/n) )

    Y = emptyList(n)

    for k in range(n//2):
        Y[k] = t[k]+W*u[k]
        Y[k+n//2] = t[k] - W * u[k]
        W = W*W_n

    return Y

#Reversed FFT
def iFFT(n,Y: list):
    def RecFFT_ReverseW(n: int, a: list,original_n: int):
        spaces = original_n-n                  #for print

        if n == 1:
            return [a[0]]

        t = RecFFT_ReverseW(n // 2, [a[i] for i in range(0, n, 2)],original_n)
        u = RecFFT_ReverseW(n // 2, [a[i] for i in range(1, n, 2)],original_n)

        e = np.e;pi = np.pi
        W = 1
        W_n = e ** (-1j * ((2 * pi) / n))  # changed: POW -1

        Y = emptyList(n)

        for k in range(n // 2):
            Y[k] = t[k] + W * u[k]
            Y[k + n // 2] = t[k] - W * u[k]
            W = W * W_n

        return Y

    v = RecFFT_ReverseW(n,Y,n)
    a = [v_i/n for v_i in v]
    return a

def mulPol(a: list, b: list):
    n = abFixLen_and_get_n(a, b) #[1,1,3] -> [1,1,3,0] , [-1,1] -> [-1,1,0,0] , n=4
    Ya = FFT(n,a)
    Yb = FFT(n,b)
    Yc = [Ya[i]*Yb[i] for i in range(n)]
    #print("\nYc =", pretty_str(Yc))
    c = iFFT(n,Yc)
    return c


# Multiplying 2 Polynoms in efficiency by FFT - O(nlogn) instead of n^2 (regular algo) :
# c(x) = a(x) * b(x)
a=[1,1,3]            #a(x)=3x^2+1x^1+1
b=[-1,1]            #b(x)=1x^1+-1

c=mulPol(a,b)       #c(x) = =3x^3+-2x^2+x^1+-1
print_polynom(c,'c')
