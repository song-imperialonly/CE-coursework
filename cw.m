%
% Readme:
% 
% These code run in the MATLAB R2025b environment. This code file is 
% a collection of separate code used to solve all questions in part B.
% 
% Thus, please do not run the entire code file directly.
% 
% If necessary, please extract the code corresponding to each question 
% and run it in conjunction with the setted macro.
% 
% Author: Shirui Song
%

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

set(0, 'DefaultAxesFontSize', 32);
set(0, 'DefaultTextFontSize', 32);

clc;
clear;
close all;

% 2
M = 1; F = 1; Lp = 1; g = 9.81;
A = [0 1 0 0;
     0 -F/M 0 0;
     0 0 0 1;
     0 F/(Lp*M) g/Lp 0];
B  = [0; 1/M; 0; -1/(Lp*M)];
P1 = [0; 1/M; 0; -1/(Lp*M)];
Cy = [1 0 0 0;
      0 0 1 0];

alpha = 1;
omega = 0.1;

S = [0 0 0;
     0 0 omega;
     0 -omega 0];
P = [P1 zeros(4,2)];
Cr = [1 0 0 0];
Q  = [0 -1 0];
I4 = eye(4); I3 = eye(3);
Aeq1 = kron(S.', I4) - kron(I3, A);  
Aeq2 = -kron(I3, B);                  
beq1 = P(:);
Aeq3 = kron(I3, Cr);                
beq3 = -Q(:);
Aeq = [Aeq1 Aeq2;
       Aeq3 zeros(3,3)];
beq = [beq1;
       beq3];
theta = Aeq \ beq;
Pi    = reshape(theta(1:12),4,3);
Gamma = theta(13:end).';  
p_des = [-0.8 -1.1 -1.4 -1.7];
K = place(A, -B, p_des);
poles_cl = eig(A + B*K);

Lgain = Gamma - K*Pi;      
ell1 = Lgain(1); ell2 = Lgain(2); ell3 = Lgain(3);

d1 = @(t) 0.5*square(2*pi/50*t);
d2 = @(t) alpha*sin(omega*t);
d3 = @(t) alpha*cos(omega*t);
u = @(t,x) (K*x + ell1*d1(t) + ell2*d2(t) + ell3*d3(t));
Tend = 60;
x0   = zeros(4,1);             
ode = @(t,x) (A*x + B*u(t,x) + P1*d1(t));
[t,X] = ode45(ode, [0 Tend], x0);
Y   = (Cy*X.').';
s   = Y(:,1);
phi = Y(:,2);
d2v = alpha*sin(omega*t);
d1v = 0.5*square(2*pi/50*t);
e = s - d2v;
uvec = zeros(size(t));
for k = 1:numel(t)
    uvec(k) = u(t(k), X(k,:).');
end

figure;
plot(t,d2v,'-','LineWidth',4); hold on;
plot(t,s,':','LineWidth',4);
plot(t,phi,'--','LineWidth',4);
grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$y(t)$','Interpreter','latex');
legend({'Ref $d_{2}(t)$','$s(t)$','$\phi(t)$'}, ...
       'Interpreter','latex','Location','northeast');

figure;
plot(t,uvec,'m-','LineWidth',4); grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$u(t)$','Interpreter','latex');
ylim([-1 1])

figure;
plot(t,e,'k-','LineWidth',4); grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$e(t)$','Interpreter','latex');
ylim([-1 1])

figure;
plot(t,d1v,'-','LineWidth',4); grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$d_{1}(t)$','Interpreter','latex');
ylim([-0.6 0.6])

% 3
M = 1; F = 1; Lp = 1; g = 9.81;
alpha = 1; omega = 0.1;
A = [0 1 0 0;
     0 -F/M 0 0;
     0 0 0 1;
     0 F/(Lp*M) g/Lp 0];
B  = [0; 1/M; 0; -1/(Lp*M)];
P1 = [0; 1/M; 0; -1/(Lp*M)];
Cy = [1 0 0 0;
      0 0 1 0];
S = [0 0 0;
     0 0 omega;
     0 -omega 0];
P = [P1 zeros(4,2)];
Cr = [1 0 0 0];
Q  = [0 -1 0];
I4 = eye(4); I3 = eye(3);
Aeq1 = kron(S.', I4) - kron(I3, A);
Aeq2 = -kron(I3, B);
beq1 = P(:);
Aeq3 = kron(I3, Cr);
beq3 = -Q(:);
Aeq = [Aeq1 Aeq2;
       Aeq3 zeros(3,3)];
beq = [beq1;
       beq3];
theta = Aeq \ beq;
Pi    = reshape(theta(1:12),4,3);
Gamma = theta(13:end).';
p_des = [-0.8 -1.1 -1.4 -1.7];
K = place(A, -B, p_des);         
Lgain = Gamma - K*Pi;             
ell1 = Lgain(1); ell2 = Lgain(2); ell3 = Lgain(3);
d1 = @(t) 0.5*square(2*pi/50*t);   
d2 = @(t) alpha*sin(omega*t);
d3 = @(t) alpha*cos(omega*t);
u = @(t,x) (K*x + ell1*d1(t) + ell2*d2(t) + ell3*d3(t));
f_nl = @(t,x) [
    x(2);
    -F*x(2) + u(t,x) + d1(t);                                        
    x(4);
    g*sin(x(3)) - (-F*x(2) + u(t,x) + d1(t))*cos(x(3))];
Tend = 60;
x0   = zeros(4,1);
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
[t,X] = ode45(f_nl, [0 Tend], x0, opts);
Y   = (Cy*X.').';
s   = Y(:,1);
phi = Y(:,2);
d2v = alpha*sin(omega*t);
d1v = 0.5*square(2*pi/50*t);
uvec = zeros(size(t));
for k = 1:numel(t)
    uvec(k) = u(t(k), X(k,:).');
end
e = s - d2v;

figure;
plot(t,d2v,'-','LineWidth',4); hold on;
plot(t,s,':','LineWidth',4);
plot(t,phi,'--','LineWidth',4);
grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$y(t)$','Interpreter','latex');
legend({'Ref $d_{2}(t)$','$s(t)$','$\phi(t)$'}, ...
       'Interpreter','latex','Location','northeast');

figure;
plot(t,uvec,'m-','LineWidth',4); grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$u(t)$','Interpreter','latex');
ylim([-1 1])

figure;
plot(t,e,'k-','LineWidth',4); grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$e(t)$','Interpreter','latex');
ylim([-1 1])

figure;
plot(t,d1v,'-','LineWidth',4); grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$d_{1}(t)$','Interpreter','latex');
ylim([-0.6 0.6])

% 4
M = 1; F = 1; Lp = 1; g = 9.81;
A = [0 1 0 0;
     0 -F/M 0 0;
     0 0 0 1;
     0 F/(Lp*M) g/Lp 0];
B  = [0; 1/M; 0; -1/(Lp*M)];
P1 = [0; 1/M; 0; -1/(Lp*M)];
Cy = [1 0 0 0;
      0 0 1 0];
alpha = 1;
omega = 1;
S = [0 0 0;
     0 0 omega;
     0 -omega 0];
P = [P1 zeros(4,2)];
Cr = [1 0 0 0];
Q  = [0 -1 0];
I4 = eye(4); I3 = eye(3);
Aeq1 = kron(S.', I4) - kron(I3, A);   
Aeq2 = -kron(I3, B);                 
beq1 = P(:);
Aeq3 = kron(I3, Cr);                 
beq3 = -Q(:);
Aeq = [Aeq1 Aeq2;
       Aeq3 zeros(3,3)];
beq = [beq1;
       beq3];
theta = Aeq \ beq;
Pi    = reshape(theta(1:12),4,3);
Gamma = theta(13:end).';  
p_des = [-0.5 -0.8 -1.1 -1.4];
K = place(A, -B, p_des);
poles_cl = eig(A + B*K);

Lgain = Gamma - K*Pi;    
ell1 = Lgain(1); ell2 = Lgain(2); ell3 = Lgain(3);
d1 = @(t) 0.5*square(2*pi/50*t);
d2 = @(t) alpha*sin(omega*t);
d3 = @(t) alpha*cos(omega*t);
u = @(t,x) (K*x + ell1*d1(t) + ell2*d2(t) + ell3*d3(t));
Tend = 20;
x0   = zeros(4,1);               
ode = @(t,x) (A*x + B*u(t,x) + P1*d1(t));
[t,X] = ode45(ode, [0 Tend], x0);
Y   = (Cy*X.').';
s   = Y(:,1);
phi = Y(:,2);
d2v = alpha*sin(omega*t);
d1v = 0.5*square(2*pi/50*t);
e = s - d2v;
uvec = zeros(size(t));
for k = 1:numel(t)
    uvec(k) = u(t(k), X(k,:).');
end

figure;
plot(t,d2v,'-','LineWidth',4); hold on;
plot(t,s,':','LineWidth',4);
plot(t,phi,'--','LineWidth',4);
grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$y(t)$','Interpreter','latex');
legend({'Ref $d_{2}(t)$','$s(t)$','$\phi(t)$'}, ...
       'Interpreter','latex','Location','southeast');

figure;
plot(t,uvec,'m-','LineWidth',4); grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$u(t)$','Interpreter','latex');

M = 1; F = 1; Lp = 1; g = 9.81;
alpha = 1; omega = 1;
A = [0 1 0 0;
     0 -F/M 0 0;
     0 0 0 1;
     0 F/(Lp*M) g/Lp 0];
B  = [0; 1/M; 0; -1/(Lp*M)];
P1 = [0; 1/M; 0; -1/(Lp*M)];
Cy = [1 0 0 0;
      0 0 1 0];
S = [0 0 0;
     0 0 omega;
     0 -omega 0];
P = [P1 zeros(4,2)];
Cr = [1 0 0 0];
Q  = [0 -1 0];
I4 = eye(4); I3 = eye(3);
Aeq1 = kron(S.', I4) - kron(I3, A);
Aeq2 = -kron(I3, B);
beq1 = P(:);
Aeq3 = kron(I3, Cr);
beq3 = -Q(:);
Aeq = [Aeq1 Aeq2;
       Aeq3 zeros(3,3)];
beq = [beq1;
       beq3];
theta = Aeq \ beq;
Pi    = reshape(theta(1:12),4,3);
Gamma = theta(13:end).';
p_des = [-0.5 -0.8 -1.1 -1.4];
K = place(A, -B, p_des);           
Lgain = Gamma - K*Pi;            
ell1 = Lgain(1); ell2 = Lgain(2); ell3 = Lgain(3);
d1 = @(t) 0.5*square(2*pi/50*t);   
d2 = @(t) alpha*sin(omega*t);
d3 = @(t) alpha*cos(omega*t);
u = @(t,x) (K*x + ell1*d1(t) + ell2*d2(t) + ell3*d3(t));
f_nl = @(t,x) [
    x(2);
    -F*x(2) + u(t,x) + d1(t);                                       
    x(4);
    g*sin(x(3)) - (-F*x(2) + u(t,x) + d1(t))*cos(x(3))];
Tend = 20;
x0   = zeros(4,1);
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
[t,X] = ode45(f_nl, [0 Tend], x0, opts);
Y   = (Cy*X.').';
s   = Y(:,1);
phi = Y(:,2);
d2v = alpha*sin(omega*t);
d1v = 0.5*square(2*pi/50*t);
uvec = zeros(size(t));
for k = 1:numel(t)
    uvec(k) = u(t(k), X(k,:).');
end
e = s - d2v;

figure;
plot(t,d2v,'-','LineWidth',4); hold on;
plot(t,s,':','LineWidth',4);
plot(t,phi,'--','LineWidth',4);
grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$y(t)$','Interpreter','latex');
legend({'Ref $d_{2}(t)$','$s(t)$','$\phi(t)$'}, ...
       'Interpreter','latex','Location','southeast');

figure;
plot(t,uvec,'m-','LineWidth',4); grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$u(t)$','Interpreter','latex');

M = 1; F = 1; Lp = 1; g = 9.81;  
A = [0 1 0 0;
     0 -F/M 0 0;
     0 0 0 1;
     0 F/(Lp*M) g/Lp 0];
B  = [0; 1/M; 0; -1/(Lp*M)];
P1 = [0; 1/M; 0; -1/(Lp*M)];
Cy = [1 0 0 0;
      0 0 1 0];
alpha = 1;
omega = 10;
S = [0 0 0;
     0 0 omega;
     0 -omega 0];
P = [P1 zeros(4,2)];
Cr = [1 0 0 0];
Q  = [0 -1 0];
I4 = eye(4); I3 = eye(3);
Aeq1 = kron(S.', I4) - kron(I3, A);  
Aeq2 = -kron(I3, B);                
beq1 = P(:);
Aeq3 = kron(I3, Cr);                 
beq3 = -Q(:);
Aeq = [Aeq1 Aeq2;
       Aeq3 zeros(3,3)];
beq = [beq1;
       beq3];
theta = Aeq \ beq;
Pi    = reshape(theta(1:12),4,3);
Gamma = theta(13:end).';   
p_des = [-0.8 -1.1 -1.4 -1.7];
K = place(A, -B, p_des);
poles_cl = eig(A + B*K);
Lgain = Gamma - K*Pi;      
ell1 = Lgain(1); ell2 = Lgain(2); ell3 = Lgain(3);
d1 = @(t) 0.5*square(2*pi/50*t);
d2 = @(t) alpha*sin(omega*t);
d3 = @(t) alpha*cos(omega*t);
u = @(t,x) (K*x + ell1*d1(t) + ell2*d2(t) + ell3*d3(t));
Tend = 20;
x0   = zeros(4,1);              
ode = @(t,x) (A*x + B*u(t,x) + P1*d1(t));
[t,X] = ode45(ode, [0 Tend], x0);
Y   = (Cy*X.').';
s   = Y(:,1);
phi = Y(:,2);
d2v = alpha*sin(omega*t);
d1v = 0.5*square(2*pi/50*t);
e = s - d2v;
uvec = zeros(size(t));
for k = 1:numel(t)
    uvec(k) = u(t(k), X(k,:).');
end

figure;
plot(t,d2v,'-','LineWidth',4); hold on;
plot(t,s,':','LineWidth',4);
plot(t,phi,'--','LineWidth',4);
grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$y(t)$','Interpreter','latex');
legend({'Ref $d_{2}(t)$','$s(t)$','$\phi(t)$'}, ...
       'Interpreter','latex','Location','northeast');

figure;
plot(t,uvec,'m-','LineWidth',4); grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$u(t)$','Interpreter','latex');

M = 1; F = 1; Lp = 1; g = 9.81;
alpha = 1; omega = 10;
A = [0 1 0 0;
     0 -F/M 0 0;
     0 0 0 1;
     0 F/(Lp*M) g/Lp 0];
B  = [0; 1/M; 0; -1/(Lp*M)];
P1 = [0; 1/M; 0; -1/(Lp*M)];
Cy = [1 0 0 0;
      0 0 1 0];
S = [0 0 0;
     0 0 omega;
     0 -omega 0];
P = [P1 zeros(4,2)];
Cr = [1 0 0 0];
Q  = [0 -1 0];
I4 = eye(4); I3 = eye(3);
Aeq1 = kron(S.', I4) - kron(I3, A);
Aeq2 = -kron(I3, B);
beq1 = P(:);
Aeq3 = kron(I3, Cr);
beq3 = -Q(:);
Aeq = [Aeq1 Aeq2;
       Aeq3 zeros(3,3)];
beq = [beq1;
       beq3];
theta = Aeq \ beq;
Pi    = reshape(theta(1:12),4,3);
Gamma = theta(13:end).';
p_des = [-0.5 -0.8 -1.1 -1.4];
K = place(A, -B, p_des);          
Lgain = Gamma - K*Pi;             
ell1 = Lgain(1); ell2 = Lgain(2); ell3 = Lgain(3);
d1 = @(t) 0.5*square(2*pi/50*t);  
d2 = @(t) alpha*sin(omega*t);
d3 = @(t) alpha*cos(omega*t);
u = @(t,x) (K*x + ell1*d1(t) + ell2*d2(t) + ell3*d3(t));
f_nl = @(t,x) [
    x(2);
    -F*x(2) + u(t,x) + d1(t);                                    
    x(4);
    g*sin(x(3)) - (-F*x(2) + u(t,x) + d1(t))*cos(x(3))];
Tend = 20;
x0   = zeros(4,1);
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
[t,X] = ode45(f_nl, [0 Tend], x0, opts);
Y   = (Cy*X.').';
s   = Y(:,1);
phi = Y(:,2);
d2v = alpha*sin(omega*t);
d1v = 0.5*square(2*pi/50*t);
uvec = zeros(size(t));
for k = 1:numel(t)
    uvec(k) = u(t(k), X(k,:).');
end
e = s - d2v;

figure;
plot(t,d2v,'-','LineWidth',4); hold on;
plot(t,s,':','LineWidth',4);
plot(t,phi,'--','LineWidth',4);
grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$y(t)$','Interpreter','latex');
legend({'Ref $d_{2}(t)$','$s(t)$','$\phi(t)$'}, ...
       'Interpreter','latex','Location','southeast');
ylim([-10 10])

figure;
plot(t,uvec,'m-','LineWidth',4); grid on;
xlabel('$t$','Interpreter','latex');
ylabel('$u(t)$','Interpreter','latex');
ylim([-16000 16000])


