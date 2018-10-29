import scipy as sp

#Datos Canal
B = 10. #mts
ss = 2. #H/V 0 si es un canal rectangula
ycini = 1. #mts
tol = 1.
yc = ycini
C0 = 1. #unidades SI, 1,49 en caso de unidades britanicas
S0 = 0.001 #Depende del slope 
n = 0.013 # Depende del material
g = 9.81
Q=30.

c=0
while tol>(10**(-1000)):
	#Calculamos Variables
	A = B*yc + ss*(yc**2)
	P = B + 2*yc*sp.sqrt(ss**2+1)
	#Q = (C0*(A**(5./3.)*(P**(-2./3.))*sp.sqrt(S0))/n)
	Tc = B + 2*ss*yc
	#Preparara Newton-Raphson
	fyc = (A**3)/(Tc)-(Q**2)/g
	dfyc = 3*A**2 - 2*ss*(A**3)/(Tc**2)
	yc2 = yc - (fyc/dfyc)

	#Calculo tolerancia
	tol = abs(yc2-yc)

	#Se actualiza yc
	yc = yc2
	c+=1

print c
print "la altura critica es: {}".format(yc)
V=Q / (B*yc + ss*(yc**2))
fr=( V ) / sp.sqrt(g*yc)
print "Froude es: {}".format(fr)
