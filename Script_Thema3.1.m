#  Topic 3.1 Different versions of the Goldbeter system, for different values of the parameters n, of the initial values and of the time t 
global vs;
global vm;
global vd;
global ks;
global k1;
global k2;
global V1;
global V2;
global V3;
global V4;
global K1;
global K2;
global K3;
global K4;
global KI;
global Km1;
global Kd;
global n;

# Case 1 - Different values in the parameters and in time t 

#Initialisation of the corresponding values to m(0),p0(0),p1(0),p2(0),pn(0)
p0=[0.1;0.25;0.25;0.25;0.25];

##Initialisation of 701 equidistant values for the variable t in space [0,70]
t=linspace(0,70,701);
#Parameter values' initialisation 
vs=1.26;
vm=0.33;
vd=1.45;
ks=0.48;
k1=1.2;
k2=1.3;
V1=2.2;
V2=1.58;
V3=4;
V4=1.5;
K1=1.5;
K2=1.3;
K3=1.4;
K4=1.4;
KI=1.7;
Km1=0.9;
Kd=0.3;
n=4;

## Solution of the system of differential equations of the Goldbeter equation using the lsode command
p=lsode("fex2",p0,t);

# Where pt symbolizes the total quantity PER protein with the formula  pt=p0+p1+p2+pn, after the solution of the system of Goldbeter differential equations 
pt=p(:,2)+p(:,3)+p(:,4)+p(:,5);

## Graph of the values of P protein (p0,p1,p2,pn,pt)concentrations as a function of time with different values for the under-study case 1
figure 1;
PLOT1=plot(t,p(:,2),'b','linewidth',4);
hold on;
PLOT2=plot(t,p(:,3),'r','linewidth',4);
PLOT3=plot(t,p(:,4),'g','linewidth',4);
PLOT4=plot(t,p(:,5),'m','linewidth',4);
PLOT5=plot(t,pt,'y','linewidth',4);
xlabel("Time");
ylabel("Concetrations ");
title ("Concetrations of P vs Time");
legend([PLOT1,PLOT2,PLOT3,PLOT4,PLOT5],'p0', 'p1', 'p2','nuclear PER protein (pn)','PER protein (pt)', 'northeast' );
hold off;

## Graph of the values m, pn, pt as a function of time 
figure (2);
PLOT6=plot(t,p(:,1),'y','linewidth',4);
hold on;
PLOT7=plot(t,p(:,5),'c','linewidth',4);
PLOT8=plot(t,pt,'d','linewidth',4);
xlabel("Time");
ylabel("Concetrations of P");
title ("The simulated Concetrations of mRNA(m), PER protein (pt), nuclear PER protein (pn)) over Time");
legend([PLOT6,PLOT7,PLOT8],'mRNA(m)','nuclear PER protein (pn)','PER protein (pt)', 'northeast' );
hold off;

##  Graph Representation of PHASE SPACE of values m, pt  
figure (3);
plot(p(:,1),pt);
xlabel("mRNA(m)");
ylabel("PER protein (pt)");
title ("Phase Spaece Diagram of mRNA(m) - PER protein (pt)");
grid on;

# Case 2 - Different values in the initial values and in time t 

#Initialisation of the corresponding values to m(0),p0(0),p1(0),p2(0),pn(0)
p0=[1.5;0.34;0.78;0.63;0.99];

##Initialisation of 101 equidistant values for the variable t in space [0,100]
t=linspace(0,100,1001);
#Parameter values' initialisation 
vs=0.76;
vm=0.65;
vd=0.95;
ks=0.38;
k1=1.9;
k2=1.3;
V1=3.2;
V2=1.58;
V3=5;
V4=2.5;
K1=1;
K2=1;
K3=1;
K4=1;
KI=1;
Km1=0.5;
Kd=0.2;
n=4;

## Solution of the system of differential equations of the Goldbeter equation using the lsode command
p=lsode("fex2",p0,t);

# Where pt symbolizes the total quantity PER protein with the formula  pt=p0+p1+p2+pn, after the solution of the system of Goldbeter differential equations 
pt=p(:,2)+p(:,3)+p(:,4)+p(:,5);

## Graph of the values of P protein (p0,p1,p2,pn,pt) concentrations as a function of time with different values for the under-study case 1
figure (4);
PLOT1=plot(t,p(:,2),'b','linewidth',4);
hold on;
PLOT2=plot(t,p(:,3),'r','linewidth',4);
PLOT3=plot(t,p(:,4),'g','linewidth',4);
PLOT4=plot(t,p(:,5),'m','linewidth',4);
PLOT5=plot(t,pt,'y','linewidth',4);
xlabel("Time");
ylabel("Concetrations ");
title ("Concetrations of P vs Time");
legend([PLOT1,PLOT2,PLOT3,PLOT4,PLOT5],'p0', 'p1', 'p2','nuclear PER protein (pn)','PER protein (pt)', 'northeast' );
hold off;

## Graph of the values m, pn, pt as a function of time 
figure (5);
PLOT6=plot(t,p(:,1),'y','linewidth',4);
hold on;
PLOT7=plot(t,p(:,5),'c','linewidth',4);
PLOT8=plot(t,pt,'d','linewidth',4);
xlabel("Time");
ylabel("Concetrations of P");
title ("The simulated Concetrations of mRNA(m), PER protein (pt), nuclear PER protein (pn)) over Time");
legend([PLOT6,PLOT7,PLOT8],'mRNA(m)','nuclear PER protein (pn)','PER protein (pt)', 'northeast' );
hold off;

##  Graph Representation of PHASE SPACE of values m, pt   
figure (6);
plot(p(:,1),pt);
xlabel("mRNA(m)");
ylabel("PER protein (pt)");
title ("Phase Spaece Diagram of mRNA(m) - PER protein (pt)");
grid on;

# Case 3 - Different values in initial values, parameters and time t 

#Initialisation of the corresponding values to m(0),p0(0),p1(0),p2(0),pn(0)
p0=[1.5;0.34;0.78;0.63;0.99];

##Initialisation of 701 equidistant values for the variable t in space [0,70]
t=linspace(0,70,701);
#Parameter values' initialisation 
vs=1.26;
vm=0.33;
vd=1.45;
ks=0.48;
k1=1.2;
k2=1.3;
V1=2.2;
V2=1.58;
V3=4;
V4=1.5;
K1=1.5;
K2=1.3;
K3=1.4;
K4=1.4;
KI=1.7;
Km1=0.9;
Kd=0.3;
n=4;

## Solution of the system of differential equations of the Goldbeter equation using the lsode command
p=lsode("fex2",p0,t);

# Where pt symbolizes the total quantity PER protein with the formula  pt=p0+p1+p2+pn, after the solution of the system of Goldbeter differential equations 
pt=p(:,2)+p(:,3)+p(:,4)+p(:,5);

## Graph of the values of P protein (p0,p1,p2,pn,pt) concentrations as a function of time with different values for the under-study case 1
figure (7);
PLOT1=plot(t,p(:,2),'b','linewidth',4);
hold on;
PLOT2=plot(t,p(:,3),'r','linewidth',4);
PLOT3=plot(t,p(:,4),'g','linewidth',4);
PLOT4=plot(t,p(:,5),'m','linewidth',4);
PLOT5=plot(t,pt,'y','linewidth',4);
xlabel("Time");
ylabel("Concetrations ");
title ("Concetrations of P vs Time");
legend([PLOT1,PLOT2,PLOT3,PLOT4,PLOT5],'p0', 'p1', 'p2','nuclear PER protein (pn)','PER protein (pt)', 'northeast' );
hold off;

## Graph of the values m, pn, pt as a function of time 
figure (8);
PLOT6=plot(t,p(:,1),'y','linewidth',4);
hold on;
PLOT7=plot(t,p(:,5),'c','linewidth',4);
PLOT8=plot(t,pt,'d','linewidth',4);
xlabel("Time");
ylabel("Concetrations of P");
title ("The simulated Concetrations of mRNA(m), PER protein (pt), nuclear PER protein (pn)) over Time");
legend([PLOT6,PLOT7,PLOT8],'mRNA(m)','nuclear PER protein (pn)','PER protein (pt)', 'northeast' );
hold off;

##  Graph Representation of PHASE SPACE of values m, pt 
figure (9);
plot(p(:,1),pt);
xlabel("mRNA(m)");
ylabel("PER protein (pt)");
title ("Phase Spaece Diagram of mRNA(m) - PER protein (pt)");
grid on;