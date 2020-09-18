#author: Luiz Duczmal duczmal@ufmg.br 2020/04/05 

escalay = 1400000
nmax = 180
Ntotal = 2500000
fracg = c( 0.1170, 0.2176, 0.5233, 0.1421 )
ng = length( fracg )
nada = 0
vertical = 1

if( vertical == 0 ) {
    f=0.25
    F = matrix(
      c( f,f,f,f,
         f,f,f,f,
         f,f,f,f,
         f,f,f,f ),
        nrow = ng,
        ncol = ng,               
        byrow = FALSE 
    )    
}

if( vertical == 1 ) {
    f = 0.25
    F = matrix(
        c( 1,1,1,f,
           1,1,1,f,    
           1,1,1,f,
           f,f,f,f),
        nrow=ng,               
        ncol=ng,               
        byrow = FALSE
    )   
}   

tt = rep( 0, nmax )
h = 1.0
h2 = h/2
N = rep( Ntotal,ng )
N = N*fracg
beta = 1.23
Dr = rep( 3.5, ng ) #10.0 #3.5
Dn = rep( 3.5, ng ) #10.0 #3.5
alpha = 0.05    #0.01 
mu = 1.00
Z = 3.7  #5.2

S = rep( 0, ng )
E = rep( 0, ng )
I = rep( 0, ng )
R = rep( 0, ng )
Ir = rep( 0, ng )
In = rep( 0, ng )

S1 = rep( 0, ng )
E1 = rep( 0, ng )
I1 = rep( 0, ng )
R1 = rep( 0, ng )
Ir1 = rep( 0, ng )
In1 = rep( 0, ng )
  
K1s = rep( 0, ng )
K1E = rep( 0, ng )
k1Ir = rep( 0, ng )
k1In = rep (0, ng )
  
K2s = rep( 0, ng )
K2E = rep( 0, ng )
k2Ir = rep( 0, ng )
k2In = rep( 0, ng )
  
K3s = rep( 0, ng )
K3E = rep( 0, ng )
k3Ir = rep( 0, ng )
k3In = rep( 0, ng )
 
K4s = rep( 0, ng )
K4E = rep( 0, ng )
k4Ir = rep(0, ng )
k4In = rep( 0, ng )
  
t0 = 1.00
tn = t0

S = rep( 0, ng )
E = rep( 0, ng )
E = fracg
S = N - fracg
I = rep( 0, ng )
R = rep( 0, ng )
R = N - ( S + E + I )
Ir = rep( 0, ng )
In = rep( 0, ng )
yStotal = rep( 0, nmax )
yEtotal = rep( 0, nmax )
yItotal = rep( 0, nmax )
yRtotal = rep( 0, nmax )
yStotal[1] = Ntotal - 1
yEtotal[1] = 1
yItotal[1] = 0
yRtotal[1] = 0
yS = matrix( 0, ng, nmax )
yE = matrix( 0, ng, nmax )
yIr = matrix( 0, ng, nmax )
yIn = matrix( 0, ng, nmax )
yI = matrix( 0, ng, nmax )
yR = matrix( 0, ng, nmax )
tt[1] = t0
yS[,1] = S
yE[,1] = E
yIr[,1] = Ir
yIn[,1] = In
yI[,1] = Ir + mu * In
yR[,1] = N - ( S + E + I )

for( i in 2:nmax ) {
  tn2 = tn + h2
  tn1 = tn + h
  
  Aux1 = ( beta ) * ( F%*%( ( ( Ir ) + mu * ( In ) )/Ntotal ) )
  k1S = h * ( -Aux1*S )
  k1E = h * ( Aux1 * S  -E/Z )
  k1Ir = h * ( alpha * E/Z - Ir/Dr )
  k1In = h * ( ( 1 - alpha ) * E/Z - In/Dn )
  
  Aux2 = ( beta ) * ( F%*%( ( ( Ir + k1Ir/2 ) + mu * ( In + k1In/2 ) )/Ntotal ) )
  k2S = h * ( -Aux2 * ( S + k1S/2 ) )
  k2E = h * ( Aux2 * ( S + k1S/2 ) - ( E + k1E/2 )/Z )
  k2Ir = h * ( alpha * (E + k1E/2 )/Z - ( Ir + k1Ir/2 )/Dr )
  k2In = h * ( ( 1 - alpha ) * ( E + k1E/2 )/Z - ( In + k1In/2 )/Dn )
  
  Aux3 = ( beta ) * ( F%*%( ( ( Ir + k2Ir/2 ) + mu * ( In + k2In/2 ) )/Ntotal ) )
  k3S = h * ( -Aux3 * ( S + k2S/2 ) )
  k3E = h * ( Aux3 * ( S + k2S/2 ) - ( E + k2E/2 )/Z )
  k3Ir = h * ( alpha * ( E + k2E/2 )/Z - ( Ir + k2Ir/2 )/Dr )
  k3In = h * ( ( 1 - alpha ) * ( E + k2E/2 )/Z - ( In + k2In/2 )/Dn )
  
  Aux4 = ( beta ) * ( F%*%( ( ( Ir + k3Ir ) + mu * ( In + k3In ) )/Ntotal ) )
  k4S = h * ( -Aux4 * ( S + k3S ) )
  k4E = h * ( Aux4 * ( S + k3S ) - ( E + k3E )/Z )
  k4Ir = h * ( alpha * ( E + k3E )/Z - ( Ir + k3Ir )/Dr )
  k4In = h * ( ( 1 - alpha ) * ( E + k3E )/Z - ( In + k3In )/Dn )
 
  S1 = S + ( k1S + 2*k2S + 2*k3S + k4S )/6
  E1 = E + ( k1E + 2*k2E + 2*k3E + k4E )/6
  Ir1 = Ir + ( k1Ir + 2*k2Ir + 2*k3Ir + k4Ir )/6
  In1 = In + ( k1In + 2*k2In + 2*k3In + k4In )/6
  
  tt[i] = tn1
  yS[,i] = S1
  yE[,i] = E1
  yIr[,i] = Ir1
  yIn[,i] = In1
  yI[,i] = Ir1 + mu*In1
  yR[,i] = N - ( S1 + E1 + Ir1 + mu*In1 )
  yStotal[i] = sum( S1 )
  yEtotal[i] = sum( E1 )
  yItotal[i] = sum( Ir1 + mu*In1 )
  yRtotal[i] = sum( N - ( S1 + E1 + Ir1 + mu*In1 ) )
  tn = tn1
  S = S1
  E = E1
  Ir = Ir1
  In = In1
  I = Ir1 + mu*In1
  R = N - ( S1 + E1 + Ir1 + mu*In1 )
}

if( ( vertical == 0 ) && ( nada == 0 ) ) {
  titulo = paste( 'Age groups: 0-9 (12%)  10-24 (22%)  25-59 (52%)  60+ (14%)',
                  '\n',
                  'Horizontal Isolation: 4-fold reduction for all age groups',
                  sep = '' )
}

if( vertical == 1 ) { 
    titulo = paste( 'Age groups: 0-9 (12%)  10-24 (22%)  25-59 (52%)  60+ (14%)',
                    '\n',
                    'Vertical Isolation: 4-fold reduction only for 60+ years age group',
                    sep = '' )
}

if( nada == 1 ) {
    titulo = paste( 'Age groups: 0-9 (12%)  10-24 (22%)  25-59 (52%)  60+ (14%)',
                    '\n',
                    'No intervention',
                    sep = '' )
}

options( scipen = 7 )
#print( yR )
#write.table( yR, file="./data/seir.csv" )
