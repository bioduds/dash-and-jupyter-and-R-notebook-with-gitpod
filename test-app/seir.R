#author: Luiz Duczmal duczmal@ufmg.br 2020/04/05 

escalay=1400000
nmax=180
Ntotal=2500000
#fracg=c(0.0565,0.0605,0.0695,0.0708,0.0773,0.0864,0.0962,0.087,0.069,0.0647,0.0647,0.0553,0.0447,0.034,0.0233,0.0176,0.0225)
#fracg=c(0.117,0.1403,0.1637,0.1832,0.1337,0.12,0.0787,0.0409,0.0225)
fracg=c(0.1170,0.2176,0.5233,0.1421)
ng=length(fracg)

nada=0
vertical=1
if(vertical==0){
 #f=1
 f=0.25
 #f=0.55
 F = matrix(   
  c(f,f,f,f,
    f,f,f,f,
    f,f,f,f,
    f,f,f,f),
   nrow=ng,               
   ncol=ng,               
   byrow = FALSE)    
}
if(vertical==1){
 f=0.25
 F = matrix( 
  c(1,1,1,f,
    1,1,1,f,
    1,1,1,f,
    f,f,f,f),
   nrow=ng,               
   ncol=ng,               
   byrow = FALSE)   
}   

tt=rep(0,nmax)
h=1.0
h2=h/2
N=rep(Ntotal,ng)
N=N*fracg
beta=1.23
Dr=rep(3.5,ng) #10.0 #3.5
Dn=rep(3.5,ng) #10.0 #3.5
alpha=0.05    #0.01 
mu=1.00
Z=3.7  #5.2

S=rep(0,ng)
E=rep(0,ng)
I=rep(0,ng)
R=rep(0,ng)
Ir=rep(0,ng)
In=rep(0,ng)

S1=rep(0,ng)
E1=rep(0,ng)
I1=rep(0,ng)
R1=rep(0,ng)
Ir1=rep(0,ng)
In1=rep(0,ng)
  
K1s=rep(0,ng)
K1E=rep(0,ng)
k1Ir=rep(0,ng)
k1In=rep(0,ng)
  
K2s=rep(0,ng)
K2E=rep(0,ng)
k2Ir=rep(0,ng)
k2In=rep(0,ng)
  
K3s=rep(0,ng)
K3E=rep(0,ng)
k3Ir=rep(0,ng)
k3In=rep(0,ng)
  
K4s=rep(0,ng)
K4E=rep(0,ng)
k4Ir=rep(0,ng)
k4In=rep(0,ng)
  
t0=1.00
tn=t0


S=rep(0,ng)
E=rep(0,ng)
E=fracg
S=N-fracg
I=rep(0,ng)
R=rep(0,ng)
R=N-(S+E+I)
Ir=rep(0,ng)
In=rep(0,ng)
yStotal=rep(0,nmax)
yEtotal=rep(0,nmax)
yItotal=rep(0,nmax)
yRtotal=rep(0,nmax)
yStotal[1]=Ntotal-1
yEtotal[1]=1
yItotal[1]=0
yRtotal[1]=0
yS=matrix(0,ng,nmax)
yE=matrix(0,ng,nmax)
yIr=matrix(0,ng,nmax)
yIn=matrix(0,ng,nmax)
yI=matrix(0,ng,nmax)
yR=matrix(0,ng,nmax)
tt[1]=t0
yS[,1]=S
yE[,1]=E
yIr[,1]=Ir
yIn[,1]=In
yI[,1]=Ir+mu*In
yR[,1]=N-(S+E+I)

for(i in 2:nmax){
  tn2=tn+h2
  tn1=tn+h
  
  Aux1=(beta)*(F%*%(((Ir)+mu*(In))/Ntotal))
  k1S=h*(   -Aux1*S              )
  k1E=h*(    Aux1*S  -E/Z        )
  k1Ir=h*(    alpha *E/Z - Ir/Dr )
  k1In=h*( (1-alpha)*E/Z - In/Dn )
  
  Aux2=(beta)*(F%*%(((Ir+k1Ir/2)+mu*(In+k1In/2))/Ntotal))
  k2S=h*( -Aux2*(S+k1S/2)                               )
  k2E=h*(  Aux2*(S+k1S/2)  -(E+k1E/2)/Z                 )
  k2Ir=h*(          alpha *(E+k1E/2)/Z - (Ir+k1Ir/2)/Dr )
  k2In=h*(       (1-alpha)*(E+k1E/2)/Z - (In+k1In/2)/Dn )
  
  Aux3=(beta)*(F%*%(((Ir+k2Ir/2)+mu*(In+k2In/2))/Ntotal))
  k3S=h*( -Aux3*(S+k2S/2)                               )
  k3E=h*(  Aux3*(S+k2S/2)  -(E+k2E/2)/Z                 )
  k3Ir=h*(           alpha*(E+k2E/2)/Z - (Ir+k2Ir/2)/Dr )
  k3In=h*(       (1-alpha)*(E+k2E/2)/Z - (In+k2In/2)/Dn )
  
  Aux4=(beta)*(F%*%(((Ir+k3Ir)+mu*(In+k3In))/Ntotal))
  k4S=h*( -Aux4*(S+k3S)                           )
  k4E=h*(  Aux4*(S+k3S)  -(E+k3E)/Z               )
  k4Ir=h*(         alpha*(E+k3E)/Z - (Ir+k3Ir)/Dr )
  k4In=h*(     (1-alpha)*(E+k3E)/Z - (In+k3In)/Dn )
 
  S1=S + (k1S +2*k2S +2*k3S +k4S )/6
  E1=E + (k1E +2*k2E +2*k3E +k4E )/6
  Ir1=Ir+(k1Ir+2*k2Ir+2*k3Ir+k4Ir)/6
  In1=In+(k1In+2*k2In+2*k3In+k4In)/6
  
  tt[i]=tn1
  yS[,i]=S1
  yE[,i]=E1
  yIr[,i]=Ir1
  yIn[,i]=In1
  yI[,i]=Ir1+mu*In1
  yR[,i]=N-(S1+E1+Ir1+mu*In1)
  yStotal[i]=sum(S1)
  yEtotal[i]=sum(E1)
  yItotal[i]=sum(Ir1+mu*In1)
  yRtotal[i]=sum(N-(S1+E1+Ir1+mu*In1))
  tn=tn1
  S=S1
  E=E1
  Ir=Ir1
  In=In1
  I=Ir1+mu*In1
  R=N-(S1+E1+Ir1+mu*In1)
}


#titulo=paste('População em isolamento: ',toString(fracg*100),'%\n Redução de contatos: ',toString(round(1/F)),'X',sep='')
if((vertical==0)&&(nada==0)){
  titulo=paste('Age groups: 0-9 (12%)  10-24 (22%)  25-59 (52%)  60+ (14%)','\n','Horizontal Isolation: 4-fold reduction for all age groups',sep='')
  #titulo=paste('Age groups: 0-9 (12%)  10-24 (22%)  25-59 (52%)  60+ (14%)','\n','Horizontal Isolation: 1.8-fold reduction for all age groups',sep='')
}
if(vertical==1)titulo=paste('Age groups: 0-9 (12%)  10-24 (22%)  25-59 (52%)  60+ (14%)','\n','Vertical Isolation: 4-fold reduction only for 60+ years age group',sep='')
if(nada==1)titulo=paste('Age groups: 0-9 (12%)  10-24 (22%)  25-59 (52%)  60+ (14%)','\n','No intervention',sep='')


#titulo=''
#Desliga o modo notação científica
options(scipen=7)
#Com escala fixa por escalay
#matplot(tt, cbind(yI[1,]+yI[2,]),type="l",col=c("black"),lty=c(1,1),main=titulo,xlab="DIAS \n a partir de 15 de março",ylab="INFECTADOS",ylim=c(0,escalay))

#matplot(tt, cbind(yI[1,],yI[2,],yI[3,],yI[4,],yItotal),type="l",col=c("blue","darkgreen","orange","red","black"),lty=c(1,1),main=titulo,xlab="DIAS \n a partir de 15 de março",ylab="INFECTADOS",ylim=c(0,escalay))
#matplot(tt, cbind(yR[1,],yR[2,],yR[3,],yR[4,],yItotal),type="l",col=c("blue","darkgreen","orange","red","black"),lty=c(1,1),main=titulo,xlab="DIAS \n a partir de 15 de março",ylab="INFECTADOS",ylim=c(0,escalay))
#matplot(tt, cbind(yR[1,],yR[2,],yR[3,],yR[4,],yR[5,],yR[6,],yR[7,],yR[8,],yR[9,],yItotal),type="l",col=c("pink","brown","red","gold","green","darkgreen","blue","violet","orange","black"),lty=c(1,1),main=titulo,xlab="DIAS \n a partir de 15 de março",ylab="INFECTADOS",ylim=c(0,escalay))
#matplot(tt, cbind(yR[1,],yR[2,],yR[3,],yR[4,],yItotal),type="l",col=c("blue","darkgreen","orange","red","black"),lty=c(1,1),main=titulo,xlab="DIAS \n a partir de 15 de março",ylab="INFECTADOS",ylim=c(0,escalay))
#matplot(tt, cbind(yR[1,],yR[2,],yR[3,],yR[4,],yItotal),type="l",col=c("blue","darkgreen","orange","red","black"),lty=c(10,3,2,5,1),lwd=3,main=titulo,xlab="DAYS \n since March 15",ylab="INFECTED PERSONS",ylim=c(0,escalay))
#abline(h=Ntotal*fracg, col=c("blue","darkgreen","orange","red"),lwd=1, lty=2)

#write.csv( tt, 'tt.csv' )
print( tt )

#Com escala automática
#matplot(tt, cbind(yy9,yy10,yy11),type="l",col=c("blue","red","black"),lty=c(1,1),main=titulo,xlab="DIAS \n a partir de 15 de março",ylab="INFECTADOS")

#legend(x=135,y=1000000,c("    0-9 accumulated","    10-24 accumulated","    25-59 accumulated","    60+ accumulated","    currently infected"),cex=0.8,col=c("blue","darkgreen","orange","red","black"),pch=c(15,15,15,15,15))