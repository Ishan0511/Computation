clc;
clear all;
syms x
z = inputdlg('Enter a single variable equation', 's');
m = cell2mat(z);
fh = str2func(['@(x)',m]);
y = fh;
yder = diff(y,x);
fhder(x) = yder;
fplot(x,y)
grid on
xlabel('x')
ylabel('f(x)')
title('Plot of the Function')
u = inputdlg({'Enter first guess value Xest:','Enter second guess value Xest2:', 'Enter Method Choice for finding the root(s): for Incremental Search Method(Type 1), for Bisection Method(Type 2), for Secant Method(Type 3), for Newton Raphson Method(Type 4) and for Regula Falsi Method(Type 5)'},'Enter Interval',1);
Xest = str2double(u{1});
Xest2 = str2double(u{2});
Choice = str2double(u{3});
fh_var = y;
fhder_var = yder;
Err = 0.0001;
Xmean = (Xest+Xest2)/2;
Xdel = Xest2-Xest;
z = fzero(fh,[Xest,Xest2]);

Xt_var = IncrementalSearch(fh,Xest,Xest2,Xdel,Err)
if Choice == 1
    Xsolution1 = Xt_var;
    msg = cell(2,1);
    msg{1} = sprintf('The Root using Incremental Search Method is :%f',Xsolution1');
    msg{2} = sprintf('The Root using fzero function is :%f',z');
    msgbox(msg)
end


Xb_var = Bisection(fh,Xest,Xest2,Xmean,Err)
if Choice == 2
    Xsolution2 = Xb_var;
    msg = cell(2,1);
    msg{1} = sprintf('The Root using Bisection Method is :%f',Xsolution2');
    msg{2} = sprintf('The Root using fzero function is :%f',z');
    msgbox(msg)
end

Xm_var = Secant(fh,Xest,Xest2,Err)
if Choice == 3
    Xsolution3 = Xm_var;
    msg = cell(2,1);
    msg{1} = sprintf('The Root using Secant Method is :%f',Xsolution3');
    msg{2} = sprintf('The Root using fzero function is :%f',z');
    msgbox(msg)   
end
    
Xs_var = NewtonRaphson(fh,fhder,Xest,Err)
if Choice == 4
    Xsolution4 = Xs_var;
    msg = cell(2,1);
    msg{1} = sprintf('The Root using Newton Raphson Method is :%f',Xsolution4');
    msg{2} = sprintf('The Root using fzero function is :%f',z');
    msgbox(msg)    
end

Xr_var = RegulaFalsi(fh,Xest,Xest2,Err)
if Choice == 5
    Xsolution5 = Xr_var;
    msg = cell(2,1);
    msg{1} = sprintf('The Root using Regula Falsi Method is :%f',Xsolution5');
    msg{2} = sprintf('The Root using fzero function is :%f',z');
    msgbox(msg)
end

function Xt = IncrementalSearch(fh,Xest,Xest2,Xdel,Err)
Xi=Xest;
for i = 1:20
    Xi = Xi + Xdel/10;
    Xt = Xi;
    if abs(fh(Xi)-fh(Xest)) < Err        
    end
    if i == 20
        fprintf('Solution was not obtained in 20 iterations,\n')
    end
    if fh(Xest)*fh(Xi)<0
        Xest = Xi-Xdel/10;
        Xdel = Xi - Xest;
        Xi=Xest;
    end
end
end
        
    
function Xb = Bisection(fh,Xest,Xest2,Xmean,Err)
while abs(fh(Xmean))>Err
    if (fh(Xmean)*fh(Xest2))<0
        Xest = Xmean;
    else
        Xest2 = Xmean;
    end
    Xmean = (Xest+Xest2)/2;
    Xb = Xmean;
end
end


function Xm = Secant(fh,Xest,Xest2,Err)
for i = 1:20
    Xi = Xest2 - (fh(Xest2)*(Xest2-Xest))/(fh(Xest2)-fh(Xest))
    if abs((Xi-Xest2)/Xest2)<Err
        Xm = Xi;
        break
    end
    Xest2 = Xi;
    if i == 20
        fprintf('Solution was not obtained in 20 iterations,\n')
        Xm = ('No answer');
    end
end
end
    
function Xs = NewtonRaphson(fh,fhder,Xest,Err)
for i = 1:20
    Xi = Xest-fh(Xest)/double(fhder(Xest))
    if abs((Xi-Xest)/Xest)<Err
        Xs = Xi;
        break
    end
    Xest = Xi;
    if i == 20
        fprintf('Solution was not obtained in 20 iterations,\n')
        Xs = ('No answer');
    end
end
end

function Xr = RegulaFalsi(fh,Xest,Xest2,Err)
for i = 1:20
    Xi = Xest2 - (fh(Xest2)*(Xest2-Xest))/(fh(Xest2)-fh(Xest))
    if abs((Xi-Xest2)/Xest2)<Err
        Xr = Xi;
        break
    end
    Xest2 = Xi;
    if i == 20
        fprintf('Solution was not obtained in 20 iterations,\n')
        Xr = ('No answer');
    end
end
end




    