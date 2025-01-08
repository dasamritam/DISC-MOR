clear;
clc;
close all;

%% set parameters
set_parameters;

%% set temperature
T0 = 30;
TM = 30;

%% set spatial domain
N = 1000;
R_find = 50;
R_use = 50;
R_POD = 10;
X = linspace(S(1),S(end),N); dx = X(2) - X(1);
T = linspace(0,1e-3,200);

%% set input
kcal = 5e4;
input1=kcal*square(5000*pi*T);
input2=-kcal*square(5000*pi*T);

u = [input1; input2];

%% construct the intial condition and steady-state solution

T_0 = zeros(1,N);
T_st = zeros(1,N);

for i = 1:N
    x = X(i);
    T_0(i) = T_i0(x,S);
    T_st(i) = T0*g_i1(x,K,H,S) + TM*g_i2(x,K,H,S);
end

plot_wall(S,100,10); 
hold on;
plot_temperature_profile(X,T_0,'k:');
plot_temperature_profile(X,T_st,'k--');
plot_labels(85,TM,S);
axis tight

%% construct Gamma
syms omega
[Gamma,sym] = construct_GAMMA();

SYM = [S(1), S(2), S(3), S(4),...
       H(1), H(2), H(3), H(4),...
       K(1), K(2), K(3),...
       D(1), D(2), D(3), omega];
   
Gamma = subs(Gamma,sym,SYM);

%% find omega_n and function U_i

F = simplify(det(Gamma));
f = matlabFunction(F);

x0 = 0; 
Omega = [0];
Phi = zeros(R_find,N);
Phin = zeros(R_find,N);

i = 2;

% create progress bar
text = ['Finding ', num2str(R_find), ' natural modes']; h = waitbar(0,text);

while length(Omega) <= R_find
    
    [Omega_, Phi_, Phin_, Index] = find_omega(f,x0,omega,Gamma,X,S,D,N);
    
    if min(abs(Omega - Omega_)) > 0.5
        Omega(i) = Omega_;
        Phi(i,:) = Phi_;
        Phin(i,:) = Phin_;
        i = i+1;
    end
    
    % update progess bar
    waitbar(length(Omega)/R_find);
    x0 = x0+1;
end

close(h) 

I1 = find(Index==1);
I2 = find(Index==2);
I3 = find(Index==3);

% delete omega = 0, and phi(0,:)
Omega(1) = [];
Phi(1,:) = [];
Phin(1,:) = [];

% sort the natural frequencies
[Omega, ind] = sort(Omega);
Phi = Phi(ind,:);

% plot first 18 natural frequencies
figure;
plot_natural_modes(Phi,Omega,X,S,3,4,'k-');
figure;
plot_natural_modes(Phin,Omega,X,S,3,4,'k-');

%% check orthogonality
R = R_use;

for i = 1:R
    for j = 1:R
        A(i,j) = trapz(X,Phi(i,:).*Phi(j,:));
        An(i,j) = trapz(X,Phin(i,:).*Phin(j,:));
    end
end

figure
subplot(1,2,1);
surf(A); axis tight; view(45,30);
title('$<\phi_n(x),\phi_k(x)>$','interpreter','latex')
subplot(1,2,2);
surf(An); axis tight; view(45,30);
title('$<\varphi_n(x),\varphi_k(x)>$','interpreter','latex')

%% fit intial condition
w_0 = T_0 - T_st;

c = zeros(R,N);
cn = zeros(R,N);
C = zeros(R,R);
Cn = zeros(R,R);

for i = 1:R
    for j = 1:R
        C(i,j) = (trapz(X,Phi(i,:).*Phi(j,:)));
        Cn(i,j) = (trapz(X,Phin(i,:).*Phin(j,:)));
    end
        B(i,1) = trapz(X,w_0.*Phi(i,:));
        Bn(i,1) = trapz(X,w_0.*Phin(i,:));
end

c = linsolve(C,B); 
cn = linsolve(Cn,Bn); 

U = (repmat(c,1,N).*Phi).'*ones(R,1);
Un = (repmat(cn,1,N).*Phin).'*ones(R,1);

figure(1);
plot_temperature_profile(X,U.' + T_st,'b-');

%% construct state space

l1 = ones(size(X)); l1(I2) = 0; l1(I3) = 0;
l2 = ones(size(X)); l2(I1) = 0; l1(I3) = 0;
l3 = ones(size(X)); l3(I2) = 0; l3(I1) = 0;
I{1} = I1; I{2} = I2; I{3} = I3; 

A = []; B = [];
for i = 1:R
    A(i,i) = -Omega(i)^2;
    B(i,1) = (D(1)^2/K(1))*trapz(X,l1.*Phin(i,:));
    B(i,2) = (D(3)^2/K(3))*trapz(X,l3.*Phin(i,:));
end

C = eye(R,R);

sys= ss(A,B,C,zeros(R,2));

%% simulate state space

[y,t,q] = lsim(sys,u,T,cn);

for i = 1:length(t)
   UR1(:,i) = (repmat(y(i,:).',1,N).*Phin).'*ones(R,1);
   UR2(:,i) = (repmat(y(i,1:3).',1,N).*Phin(1:3,:)).'*ones(3,1);
end

disp('press any key to start the simulation video')
pause;
figure(1);
for i = 1:length(T)
    cla
    plot_wall(S,100,10); 
    hold on;
    plot_temperature_profile(X,T_0,'k:');
    plot_temperature_profile(X,T_st,'k--');
    plot_temperature_profile(X,T_st+UR1(:,i).','b-');
    plot_temperature_profile(X,T_st+UR2(:,i).','r-');
    text = ['time = ',num2str(1000*T(i),'%2.3f'),' (ms)','; ','$q_1$ = ',num2str(input1(i)/1000),' (kcal)','; ','$q_2$ = ',num2str(input2(i)/1000),' (kcal)'];
    title(text,'interpreter','latex');
    text = ['FO order: ', num2str(R), '$\;$(blue); $\;$', 'FO order: ', num2str(3), '$\;$(red)'];
    xlabel(text,'interpreter','latex')
    plot_labels(85,TM,S);
    axis([S(1) S(end) 10 100]);
    pause(1/60);
end


%% POT

W = UR1*UR1.';

[U,E,V] = svd(W);

Phi_POD = (1/(sqrt(dx))*U(:,1:R_POD).');

A = [];
for i = 1:R_POD
    for j = 1:R_POD;
        A(i,j) = trapz(X,Phi_POD(i,:).*Phi_POD(j,:));
    end
end

figure(5)
subplot(1,2,1)
surf(A); axis([0 R_POD 0 R_POD -0.1 1.1]); axis tight; view(45,30);
title('$<\phi_{pod,n}(x),\phi_{pod,k}(x)>$','interpreter','latex')
subplot(1,2,2)
plot_wall(S,2,-2); hold on; 
plot(X,Phi_POD);
title('Natural modes for the POD','interpreter','latex')

%%

for i = 1:R_POD
    for j = 1:R_POD
        Cp(i,j) = (trapz(X,Phi_POD(i,:).*Phi_POD(j,:)));
    end
    
    Bp(i,1) = trapz(X,w_0.*Phi_POD(i,:));
end

cp = linsolve(Cp,Bp); 

Phi_POD_ddot=[]; 
for i = 1:R_POD
    Phi_POD_ddot=[Phi_POD_ddot; gradient(gradient(Phi_POD(i,:),dx),dx)];
end;

A = []; B = [];
for i = 1:R_POD
    for j = 1:R_POD
        A(i,j) = (D(1)^2)*trapz(X(I1),Phi_POD_ddot(i,I1).*Phi_POD(j,I1)) + (D(2)^2)*trapz(X(I2),Phi_POD_ddot(i,I2).*Phi_POD(j,I2)) + (D(3)^2)*trapz(X(I3),Phi_POD_ddot(i,I3).*Phi_POD(j,I3));
        %A(i,j) = (D(1)^2)*trapz(X,Phi_POD_ddot(i,:).*Phi_POD(j,:)) + (D(2)^2)*trapz(X,Phi_POD_ddot(i,:).*Phi_POD(j,:)) + (D(3)^2)*trapz(X,Phi_POD_ddot(i,:).*Phi_POD(j,:));
    end
    
    B(i,1) = (D(1)^2/K(1))*trapz(X,l1.*Phi_POD(i,:));
    B(i,2) = (D(3)^2/K(3))*trapz(X,l3.*Phi_POD(i,:));
end

C = eye(R_POD,R_POD);
D = zeros(R_POD,2);

sys = ss(A.',B,C,D);

[y,t,q] = lsim(sys,u,T,cp);

for i = 1:length(t)
   UPOD(:,i) = (repmat(y(i,:).',1,N).*Phi_POD).'*ones(R_POD,1);
end

disp('press any key to start the simulation video for POD')
pause;
figure(1);
for i = 1:length(T)
    cla
    plot_wall(S,100,10); 
    hold on;
    plot_temperature_profile(X,T_0,'k:');
    plot_temperature_profile(X,T_st,'k--');
    plot_temperature_profile(X,T_st+UR1(:,i).','b-');
    plot_temperature_profile(X,T_st+UR2(:,i).','g-');
    plot_temperature_profile(X,T_st+UPOD(:,i).','r-');
    text = ['time = ',num2str(1000*T(i),'%2.3f'),' (ms)','; ','$q_1$ = ',num2str(input1(i)/1000),' (kcal)','; ','$q_2$ = ',num2str(input2(i)/1000),' (kcal)'];
    title(text,'interpreter','latex');
    text = ['FO order: ', num2str(R), '$\;$(blue); $\;$', 'POD order: ', num2str(R_POD), '$\;$(red); $\;$', 'FO order: ', num2str(3), '$\;$(green)'];
    xlabel(text,'interpreter','latex')
    plot_labels(85,TM,S);
    axis([S(1) S(end) 10 100]);
    pause(1/60);
end


