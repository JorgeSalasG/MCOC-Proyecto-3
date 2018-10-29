import scipy as sp

#Datos Canal
B = 10. #mts
ss = 2. #H/V 0 si es un canal rectangula
yini = 1. #mts
tol = 1.
yn = yini
C0 = 1 #unidades SI, 1,49 en caso de unidades britanicas
S0 = 0.001 #Depende del slope 
n = 0.013 # Depende del material

while tol>(10**(-6)):
	#Calculamos Variables
	A = B*yn + ss*(yn**2)
	P = B + 2*yn*sp.sqrt(ss**2+1)
	dAdyn = B + 2*ss*yn
	dPdyn = 2*sp.sqrt(ss**2+1)
	Q = 30. #(C0*(A**(5/3)*(P**(-2/3))*sp.sqrt(S0))/n)

	#Preparara Newton-Raphson
	fyn = (A**(5/3))*(P**(-2/3))-(n*Q)/(C0*sp.sqrt(S0))
	dfyn = (5/3)*(A**(2/3))*(P**(-2/3))*dAdyn - (2/3)*(A**(5/3))*(P**(-5/3))*dPdyn
	yn2 = yn - (fyn/dfyn)

	#Calculo tolerancia
	tol = sp.absolute(yn2-yn)

	#Se actualiza yn
	yn = yn2


print "la altura normal es: {}".format(yn)
