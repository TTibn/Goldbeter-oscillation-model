## Arithmitiki epilisi toy sistimatos diaforikwn eksiswsewn tis Eksiswshs Goldbeter

##Dilwnoume global kapoies parametrous wste na allaksoume tis times toys xwris na xreiazetai na epembainoume stin function 

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

# Case 1
#Arxikopoiisi twn timwn pou antistoixoun sta m(0),p0(0),p1(0),p2(0),pn(0)
p0=[0.1;0.25;0.25;0.25;0.25];

##Arxikopoiisi 501 isapexouswn timwn gia ti metavliti t sto diastima [0,50]
t=linspace(0,50,501);
#Arxikopoiisi timwn parametrwn 
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

## Epilisi sistimatos diaforikwn eksiswsewn tis eksiswshs Goldbeter me xrisi tis entolis lsode
p=lsode("fex2",p0,t);

# Opou pt symbolizei tin oliki posotita PER protein me typo pt=p0+p1+p2+pn , meta apo ti lysi tou sistimatos diaforikwn eksiswsewn Goldbeter
pt=p(:,2)+p(:,3)+p(:,4)+p(:,5);

## Grafiki Parastasi twn timwn twn sygkentrwsewn p (p0,p1,p2,pn,pt) ws sunartisi tou xronou me diaforetikes times  gia tin ypo meleti case 1 
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

## Grafiki Parastasi twn timwn m,pn,pt ws sunartisi tou xronou 
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

## Trisdiastati Grafiki Parastasi twn timwn m,pn,pt  
figure (3);
plot(p(:,1),p(:,5),pt);
xlabel("mRNA(m)");
ylabel("nuclear PER protein (pn)");
zlabel("PER protein (pt)");
title ("3D Diagram of mRNA(m),nuclear PER protein (pn),PER protein (pt)");
grid on;

# Case 2
#Arxikopoiisi twn timwn pou antistoixoun sta m(0),p0(0),p1(0),p2(0),pn(0)
p0=[1.9;0.8;0.8;0.8;0.8];

##Arxikopoiisi 501 isapexouswn timwn gia ti metavliti t sto diastima [0,50]
t=linspace(0,50,501);
#Arxikopoiisi timwn parametrwn
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

## Epilisi sistimatos diaforikwn eksiswsewn tis eksiswshs Goldbeter me xrisi tis entolis lsode
p=lsode("fex2",p0,t);

# ���� pt ���������� ��� ����� �������� PER protein �� ���� pt=p0+p1+p2+pn , ���� ��� �� ���� ��� ��� ���������� ����. ��������� Goldbeter
pt=p(:,2)+p(:,3)+p(:,4)+p(:,5);

## Grafiki Parastasi twn timwn twn sygkentrwsewn p (p0,p1,p2,pn,pt) ws sunartisi tou xronou me diaforetikes times  gia tin ypo meleti case 2
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

## Grafiki Parastasi twn timwn m,pn,pt ws sunartisi tou xronou 
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

## Trisdiastati Grafiki Parastasi twn timwn m,pn,pt  
figure (6);
plot(p(:,1),p(:,5),pt);
xlabel("mRNA(m)");
ylabel("nuclear PER protein (pn)");
zlabel("PER protein (pt)");
title ("3D Diagram of mRNA(m),nuclear PER protein (pn),PER protein (pt)");
grid on;