clc;
clear all;
z = [0.2 ; 0.35 ; 0.45];
L = [825 ; 825 ; 1825 ; 450];
V = [0 ; 1375 ; 1375 ; 1375];
F = 1000;
D = 550;
DM3 = F*z(1,1);
DP3 = F*z(2,1);
DO3 = F*z(3,1);
D1 = [0 ; 0 ; DM3 ; 0];
D2 = [0 ; 0 ; DP3 ; 0];
D3 = [0 ; 0 ; DO3 ; 0];
x = z;
P = 29.39; %in psia
lC4 = ones(4,1);
lC5 = ones(4,1);
lC8 = ones(4,1);
%% Bubble Point Calculations

K1=@(T) exp((-1280557/(T^2))+7.94986-(0.96455*log(29.39)));
K2=@(T) exp((-1524891/(T^2))+7.33129-(0.89143*log(29.39)));
K3=@(T) exp((-7646.81641/(T))+12.48457-(0.73152*log(29.39)));
TC4=@(KC4) (1280557/(7.94986-(0.96455*log(29.39))-log(KC4)))^0.5;
TC5=@(KC5) (1524891/(7.33129-(0.89143*log(29.39))-log(KC5)))^0.5;
TC8=@(KC8) (7646.81641/(12.48457-(0.73152*log(29.39))-log(KC8)));


%% Finding Temperature by guessing K using De Priester
% We know that zc4<zc5<zc8 therefore we will guess KC8 such that KC4>KC5>1>KC8 
KC8=0.2;
Err = inf;
I = 0;  
while Err>0.01
    Temp = TC8(KC8);
    Kxsum = (K1(Temp)*z(1))+(K2(Temp)*z(2))+(K3(Temp)*z(3));
    Err=abs(Kxsum-1);
    I=I+1;
    KC8=KC8/Kxsum;
end
%% Matrix for n-C4
AC4 = zeros(4,4);
for j = 1:4 
    if j==1
        AC4(j,1) = 1+D/L(1);
        AC4(j,2) = -(K1(Temp)*V(2))/L(2);
             elseif j==4
        AC4(j,j-1) = -1;
        AC4(j,j) = 1+V(4)*K1(Temp)/L(4);
            else
        AC4(j,j-1) = -1;
        AC4(j,j) = 1+V(j)*K1(Temp)/L(j);
        AC4(j,j+1) = -(K1(Temp)*V(j+1))/L(j+1);
    end
            lC4 = AC4\D1;
end
%% Matrix for n-C5
AC5 = zeros(4,4);
for j = 1:4 
    if j==1
        AC5(j,1) = 1+D/L(1);
        AC5(j,2) = -(K2(Temp)*V(2))/L(2);
             elseif j==4
        AC5(j,j-1) = -1;
        AC5(j,j) = 1+V(4)*K2(Temp)/L(4);
            else
        AC5(j,j-1) = -1;
        AC5(j,j) = 1+V(j)*K2(Temp)/L(j);
        AC5(j,j+1) = -(K2(Temp)*V(j+1))/L(j+1);
    end
            lC5 = AC5\D2;
end

%% Matrix for n-C8
AC8 = zeros(4,4);
for j = 1:4 
    if j==1
        AC8(j,1) = 1+D/L(1);
        AC8(j,2) = -(K3(Temp)*V(2))/L(2);
             elseif j==4
        AC8(j,j-1) = -1;
        AC8(j,j) = 1+V(4)*K3(Temp)/L(4);
            else
        AC8(j,j-1) = -1;
        AC8(j,j) = 1+V(j)*K3(Temp)/L(j);
        AC8(j,j+1) = -(K3(Temp)*V(j+1))/L(j+1);
    end
            lC8 = AC8\D3;
end
%% Liquid Mole fraction for each stage
xC4=zeros(4,1);
xC5=zeros(4,1);
xC8=zeros(4,1);
for j=1:4
    xC4(j)=lC4(j)/(lC4(j)+lC5(j)+lC8(j));
    xC5(j)=lC5(j)/(lC4(j)+lC5(j)+lC8(j));
    xC8(j)=lC8(j)/(lC4(j)+lC5(j)+lC8(j));
end
%% Vapor Mole fraction and Temperature for each stage
yC4=zeros(4,1);
yC5=zeros(4,1);
yC8=zeros(4,1);
Tnew=zeros(4,1);
ITER = 0 ; 
for j=1:4
Error = inf;
Temp1 = Temp;
while Error>0.01
    Kxsum1=(K1(Temp1)*xC4(j))+(K2(Temp1)*xC5(j))+(K3(Temp1)*xC8(j));
    Error=abs(Kxsum1-1);
    ITER = ITER+1;
    KC8new=K3(Temp1)/Kxsum1;
    Temp1=TC8(KC8new);
end
yC4(j)=K1(Temp1)*xC4(j);
yC5(j)=K2(Temp1)*xC5(j);
yC8(j)=K3(Temp1)*xC8(j);
Tnew(j)=Temp1;
end
%% Temperature Convergence
MatTemp=[Temp;Temp;Temp;Temp];
I1=0;
while abs((Tnew(1)-MatTemp(1)))>0.0001 & abs((Tnew(2)-MatTemp(2)))>0.0001 & abs((Tnew(3)-MatTemp(3)))>0.0001 & abs((Tnew(4)-MatTemp(4)))>0.0001
    MatTemp(1)=Tnew(1);
    MatTemp(2)=Tnew(2);
    MatTemp(3)=Tnew(3);
    MatTemp(4)=Tnew(4);
    %% Matrix for n-C4
AC4 = zeros(4,4);
for j = 1:4 
    if j==1
        AC4(j,1) = 1+D/L(1);
        AC4(j,2) = -(K1(MatTemp(j+1))*V(2))/L(2);
             elseif j==4
        AC4(j,j-1) = -1;
        AC4(j,j) = 1+V(4)*K1(MatTemp(j))/L(4);
            else
        AC4(j,j-1) = -1;
        AC4(j,j) = 1+V(j)*K1(MatTemp(j))/L(j);
        AC4(j,j+1) = -(K1(MatTemp(j+1))*V(j+1))/L(j+1);
    end
            lC4 = AC4\D1;
end
%% Matrix for n-C5
AC5 = zeros(4,4);
for j = 1:4 
    if j==1
        AC5(j,1) = 1+D/L(1);
        AC5(j,2) = -(K2(MatTemp(j+1))*V(2))/L(2);
             elseif j==4
        AC5(j,j-1) = -1;
        AC5(j,j) = 1+V(4)*K2(MatTemp(j))/L(4);
            else
        AC5(j,j-1) = -1;
        AC5(j,j) = 1+V(j)*K2(MatTemp(j))/L(j);
        AC5(j,j+1) = -(K2(MatTemp(j+1))*V(j+1))/L(j+1);
    end
            lC5 = AC5\D2;
end

%% Matrix for n-C8
AC8 = zeros(4,4);
for j = 1:4 
    if j==1
        AC8(j,1) = 1+D/L(1);
        AC8(j,2) = -(K3(MatTemp(j+1))*V(2))/L(2);
             elseif j==4
        AC8(j,j-1) = -1;
        AC8(j,j) = 1+V(4)*K3(MatTemp(j))/L(4);
            else
        AC8(j,j-1) = -1;
        AC8(j,j) = 1+V(j)*K3(MatTemp(j))/L(j);
        AC8(j,j+1) = -(K3(MatTemp(j+1))*V(j+1))/L(j+1);
    end
            lC8 = AC8\D3;
end
%% Liquid Mole fraction for each stage
xC4=zeros(4,1);
xC5=zeros(4,1);
xC8=zeros(4,1);
for j=1:4
    xC4(j)=lC4(j)/(lC4(j)+lC5(j)+lC8(j));
    xC5(j)=lC5(j)/(lC4(j)+lC5(j)+lC8(j));
    xC8(j)=lC8(j)/(lC4(j)+lC5(j)+lC8(j));
end
    for j=1:4
Error = inf;
Temp1 = Temp;
while Error>0.01
    Kxsum1=(K1(Temp1)*xC4(j))+(K2(Temp1)*xC5(j))+(K3(Temp1)*xC8(j));
    Error=abs(Kxsum1-1);
    ITER = ITER+1;
    KC8new=K3(Temp1)/Kxsum1;
    Temp1=TC8(KC8new);
end
yC4(j)=K1(Temp1)*xC4(j);
yC5(j)=K2(Temp1)*xC5(j);
yC8(j)=K3(Temp1)*xC8(j);
Tnew(j)=Temp1;
    end
    I1=I1+1;
    Tc = (Tnew - 491.67)*5/9;
end
%% Table
Stages={'Condenser (j=1)';'Stage (j=2)';'Feed Stage (j=3)';'Reboiler (j=4)'};
xnC4=xC4;
xnC5=xC5;
xnC8=xC8;
ynC4=yC4;
ynC5=yC5;
ynC8=yC8;
Temperature=Tc;
Table=table(Stages,xnC4,xnC5,xnC8,ynC4,ynC5,ynC8,Temperature);
Calculations = rows2vars(Table)

