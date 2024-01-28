#%% Celle nr 1 
%reset -f

a = 1
print(a)

#%% Celle nr 2, feil i denne går fint for celle 1
%reset -f    

b = 2  # ufullstendig kode



#%% Presentasjon av beregningsresultat
%reset -f
import numpy as np

# legge til et tall til en variabel
a = 314.159265
print(a)

# lage strengen med + operator
str_a1 = 'a1 = ' + str(a)
print(str_a1)

# bygge opp strengen inne i print
print('a2 =',a) 

# bruke f-strenger
str_a3 = f'a3 = {a}'
print(str_a3)

# bruke f-strenger, styring av format
str_a4 = f'a4 = {a:1.2f}'
print(str_a4)


print('Element 4 er: ',str_a4[3]) 

#%% Presentasjon av beregninger, del 2
%reset -f
import numpy as np

a = 6
a += 5
a -= 3

# Modulo-divisjon
b = a%3
print('b =', b)

# liggende vektor, MÅ bruke komma ,
c = [1, 2, 3, 4]
print('c = ',c)

# transponert, d blir staaande vektor
d = np.transpose(c)
print('d =',d)




#%% Indeksering, første indeks og initialverdi
%reset -f
import numpy as np

x = np.ones(5)  
x[0] = 3
for k in range(1,5):
   x[k] = x[k-1] + 1

print(x)   

#%% Indentering, vil ikke fungere
%reset -f
import numpy as np

x = np.ones(5)  
   x[0] = 3
for k in range(1,5):
x[k] = x[k-1] + 1

print(x)   

#%% Datatyper og kommentarer
%reset -f

"""
Lange 
kommentarer 
kan 
lages slik
"""

# Datatyper
c = 3         # c blir int
d = 6.4       # d blir float
e = 5 + 3j    # e blir complex
f = int(d)    # f blir int
print('c =',type(c)) 
print('d =',type(d)) 
print('e =',type(e)) 
print('f =',type(f)) 

# Variable som starter med underscore
_a = 3
print('_a =',_a)

# Python bruker kun "j" eller "1j" som 
# kompleks variabel. 
# Kan ikke bruke "i" eller "1i" som Matlab
i1 = 5 + 3j
i2 = 5 + 3*1j

g = 5
h = 3
i5 = complex(g,h)
print(i1,i2,i5)

#%% For-løkker
%reset -f
import numpy as np

# metode 1, dynamisk utvidelse ved
# .append-funksjonen hvor A blir en liste
A = [] 
for k in range(0,10):
    A.append(k)

print('A =', A, type(A)) 

# metode 2, preallokering med np.zeros
# B blir numpy.ndarray
B = np.zeros(10)  
for k in range(0,10):
    B[k] = k

print('B =', B, type(B)) 

# metode 3, preallokering med np.ones
# hvor for kort preallokering 
# gir feilmelding
C = np.ones(10)  
for k in range(0,10):
    C[k] = k

print("C =", C, type(C))


#%% Rekursive beregninger i for-løkker
%reset -f
import numpy as np

u = np.arange(1,11)

# versjon 1)
y1 = np.sum(u)
print('y1 = ',y1)


# versjon 2)
M = len(u)
y2 = u[0]    # initialisering
for k in range(1,M):
    y2 = y2 + u[k]

print('y2 = ',y2)


# versjon 3)
y3 = np.zeros(M)
y3[0] = u[0] # initialisering
for k in range(1,M):
    y3[k] = y3[k-1] + u[k]
print('y3 = ',y3)


#%% If-setninger og logiske betingelser
%reset -f

a = input('skriv inn tallet "a" mellom 0 og 10: ')
a = float(a)

if a >= 5:
    print('  a >= 5')    
elif a > 0 and a < 5:
    print('  a > 0 og a < 5')    
else:
    print('  a <= 0')    


b = input('skriv inn tallet "b" mellom 0 og 10: ')
b = float(b)

if a == 1 or b == 9:
    print('b == 1 eller b == 9')
elif b != 1 and b != 9:
    print('  b != 1 og b != 9');






#%% List/array/vektor/matrise
%reset -f
import numpy as np

# en liste
C1 = [1, 2, 3]    # C1 = [1 2 3] gir feil

# np.arange(start,slutt,step)
C2 = np.arange(4,10,2)

"""
np.linspace(start,slutt,#verdier)
Tredje parameter er valgfri. 
Hvis utelatt så blir 
50 verdier lagt inn.
"""
C3 = np.linspace(2,8,4)

# Initialiser en 2D liste
C4 = [[1, 2, 3], [4, 5, 6]]

# Initialiser en 2D numpy.ndarray
C5 = np.array([[1, 2, 3], [4, 5, 6]])

print('C1 = ',C1,'\n',type(C1))
print('C2 = ',C2,'\n',type(C2))
print('C3 = ',C3,'\n',type(C3))
print('C4 = ',C4,'\n',type(C4))
print('C5 = ',type(C5))
print(C5)


# Skriv ut elementet med tallet 5 
# fra C4 eller C5
print('C4[1][1] = ', C4[1][1]) # rekke 2, kolonne 2

# Skriv ut rekke 1 og 2 
# fra kolonne 1 i C4.
# Skriver ut tallene 1 og 4
for i in range(0,2):
    print('C4[',i,'][0] = ',C4[i][0])





#%% Uttrekk av enkeltelement fra en liste/vektor
%reset -f
import numpy as np

# 3 alternative måter å lage liste/array
x = np.arange(1,11,1)
x = np.linspace(1,10,10)
x = list(range(1,11))
print('x =',x,'\n')

print('Første element')
print('  x[0] =',x[0],'\n')

print('Nest siste element')
print(' x[-2] =',x[-2])
print(' x[8]  =',x[8],'\n')

print('Siste element')
print(' x[-1] =',x[-1])
print(' x[9]  =',x[9],'\n')

print('Dersom x stadig blir storre,\n'
      'bruker k som argument.')
k = len(x)
print('k =',k,'\n')

print('Nest siste element')
print('x[k-2] =',x[k-2],'\n')

print('Siste element')
print('x[k-1] =',x[k-1])


#%% Slicing, uttrekk av flere elementer fra en liste/vektor
%reset -f
import numpy as np

x = np.arange(1,11,1)
print('x =',x,'\n')

print('Første tre element')
print('    x[0:3] =',x[0:3])
print('     x[:3] =',x[:3],'\n')

print('Andre og tredje element')
print('    x[1:3] =',x[1:3],'\n')

print('Tredje element')
print('    x[2:3] =',x[2:3],'\n')

print('Tredje siste og nest siste')
print('  x[-3:-1] =',x[-3:-1])
print('  x[7:9]   =',x[7:9],'\n')

print('De to siste elementene')
print('    x[-2:] =',x[-2:])
print('    x[8:]  =',x[8:],'\n')

k = 6
print('Spesifiserer vilkålig k =',k)

print('Plukker ut elementer')
print('x[k-1:]    =',x[k-1:])
print('x[k-1:k+1] =',x[k-1:k+1],'\n')


print('Indeksere utover lengden til x')
print('gir IKKE feil ved slicing!!!!!')
k = len(x)
print('Antall elementer i x, k =',k)
print('x[k-1:k+1] =',x[k-1:k+1])




#%% Stepping i en liste
%reset -f
import numpy as np

x = np.arange(1,11,1)
print('x =',x,'\n')
# x[start:stop:step]

print('Fra første element,') 
print('til, men ikke med, element Y,')
print('for hvert andre element')
print('    x[::2] =',x[::2])
print('   x[0::2] =',x[0::2])
print('  x[:-1:2] =',x[:-1:2])
print(' x[0:-1:2] =',x[0:-1:2])
print(' x[0:-2:2] =',x[0:-2:2],'\n')

print('Fra fjerde siste, til siste,')
print('for hvert andre element')
print('  x[-4::2] =',x[-4::2],'\n')

print('Fra indeks 5 til indeks 10,')
print('for hvert andre element')
print(' x[5:10:2] =',x[5:10:2],'\n')

print('Dersom x stadig blir storre,') 
print('bruker k som argument')
k = len(x)
print('Antall elementer i x, k =',k,'\n')

print('Tredje siste og siste element')
print('  x[k-3:k:2] =',x[k-3:k:2],'\n')

print('Gir IKKE feilmelding')
print('ved for høy indeks')
print('x[k-3:k+2:2] =',x[k-3:k+2:2],'\n')

print('Indekserer i ring')
print('x[k-9:k:2]  =',x[k-9:k:2])
print('x[k-10:k:2] =',x[k-10:k:2])
print('x[k-12:k:2] =',x[k-12:k:2])

#%% Plukke ut vilkårlige element i en liste
%reset -f

x = list(range(1,11))
print('x =',x,'\n')

element = [3,7,9]
x1 = [x[index] for index in element]
print('      x[3,7,9] = ',x1)

element = [-7,-3,-1]
x2 = [x[index] for index in element]
print('   x[-7,-3,-1] = ',x2)

# dersom x stadig blir storre, 
# bruker k som argument
k = len(x)

element = [k-7,k-3,k-1]
x_k = [x[index] for index in element]
print('x[k-7,k-3,k-1] = ',x_k)


#%% Plotting
%reset -f
import numpy as np
import matplotlib.pyplot as plt

Y1 = 2
Y2 = 0.5
t = np.linspace(0,2*np.pi,100)
y1 = Y1*np.sin(t)
y2 = Y2*np.cos(t)
plt.plot(t,y1,'b')
plt.plot(t,y2,'r:')
plt.plot(t,y1+y2,'g--')
plt.title('Største tid =' + f'{max(t):1.2f}') 
plt.ylabel('y-verdi')
plt.legend(['y1, ampl. Y1 = ' + str(Y1),
            'y2, ampl. Y2 = ' + str(Y2),
            'y1+y2, bare tekst'])
plt.text(1,-1,'Tekst')
plt.xlabel('tiden t [s]')
plt.show() 


#%% Egne funksjoner
%reset -f

def myfunc(arg1,arg2):
    # ** er eksponential-operator.
    # Bruker lokale variable.
    return arg1**arg2 

x = 2
y = 3
z = myfunc(x,y)
print('z =',z)  # skriver ut tallet 8 = 2^3


#%% Egne funksjoner
%reset -f

def myappendfunc(arg1,arg2):
    arg1.append(arg2)


x = [1,2,3,4]
y = 5
myappendfunc(x,y) # ingen returvariabel
print('x =',x) # skiver ut listen x med 5 elemeter.






#%%








#%%

x= [1, 2, 3, 4, 5, 6, 7,8,9,10]
k = len(x)
#print( [ x[index] for index in [4,6] ])
y=([x[index] for index in [k-3,k-1]])
     
print(y)

