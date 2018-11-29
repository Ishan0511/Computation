clc;
clear all;
T0 = 40;
T5 = 200;
A = [2.04 -1 0 0; -1 2.04 -1 0; 0 -1 2.04 -1; 0 0 -1 2.04];
B = [40.8;0.8;0.8;200.8];
Texact = A\B;
%% Square Matrix
[M,N] = size(A);
if M~=N
    error('A is not a square matrix');
end
%% Check if it is diagonally dominant
for f = 1:M
    row = abs(A(f,:));
    d = sum(row) - row(f);
    if row(f) <= d
        error('[A] is not diagonally dominant');
    end
end
%% Jacobi Iterations
tol = 1e-4;

% Initial guess
Tjacobi = zeros(M,1);

ITER = 0;
Err = inf;
while Err > tol
    Told = Tjacobi;
    for i = 1:M
        Row = 0; 
        for j = 1:M
            if j~=i
            Row = Row + A(i,j)*Told(j);   %to calculate Aij*Xj for every row
            end
                    end
        Tjacobi(i) = (B(i)-Row)/A(i,i);        
    end
    ITER = ITER + 1;
    Err = abs(Tjacobi-Told);
end

%% Gauss Siedel Iterations
tol = 1e-4;


% Initial guess
Tgs = zeros(M,1);

ITERATIONS = 0;
Err = inf;
while Err > tol
    T_old = Tgs;
    for i = 1:M
        Row = 0;
        
            for j = 1:i-1
                Row = Row + A(i,j)*Tgs(j);   %to calculate Aij*Xj for every row
            end
            for j = i+1:M
                Row = Row + A(i,j)*T_old(j);
            end
       
        
        Tgs(i) = (B(i)-Row)/A(i,i);
    end
    ITERATIONS = ITERATIONS + 1;
    Err = abs(Tgs-T_old);
end

%% Gauss Elimination
%using A(1,1) as a pivot element
f = zeros(M,1);
for i = 1:M
    f(i+1:M) = A(i+1:M,i)/A(i,i);
    for j = i+1:M
     A(j,i:M) = A(j,i:M)- f(j)*A(i,i:M);  
    end
    B(i+1:M)=B(i+1:M)-B(i)*f(i+1:M);
end
   AB = [A,B]; 

% Back- Substituition
Tge = zeros(M,1);
for i = M:-1:1
    Tge(i) = (AB(i,end) - AB(i,i+1:N)*Tge(i+1:N))/(AB(i,i));
end

%% LU Factorization
% Upper & Lower Triangle matrix [L][U] = [A]
[L,U] = lu (A);

% defining Vector Z = [U]*T      % [L]*Z = B
% calculating Z
Z = L\B;
% calculating T
Tlu = U\Z;

%% plotting the solutions
Nodes = [1 2 3 4];
plot(Nodes,Tjacobi,'x')
hold on
plot(Nodes,Tgs,'o')
plot(Nodes,Tge,'s')
plot(Nodes,Tlu,'*')
plot(Nodes,Texact,'-v')
hold off
xlabel('Nodes')
ylabel('Temperatures')
title('Plot of Temperatures at different Nodes')
legend('Tjacobi','Tgs','Tge','Tlu','Texact')

%% printing results in a table
filename = 'Temperature_Table.xlsx';
Table = {'Gauss Elimination','LU Factorization','Jacobi Iteration','Gauss Siedel Iteration','Exact Solutions'; Tge(1),Tlu(1),Tjacobi(1),Tgs(1),Texact(1);Tge(2),Tlu(2),Tjacobi(2),Tgs(2),Texact(2);Tge(3),Tlu(3),Tjacobi(3),Tgs(3),Texact(3);Tge(4),Tlu(4),Tjacobi(4),Tgs(4),Texact(4)} 
sheet = 1;
xlRange = 'B2';
xlswrite(filename,Table,sheet,xlRange)


