function [Omega, Phi, Phin, index] = find_omega(f,x0,omega,Gamma,X,S,D,N)

options = optimset('TolX',1e-6);
Omega = fzero(f,x0,options);
G1 = double(subs(Gamma,omega,Omega));

[U,E,V] = svd(G1);

sol = V(:,end);

A = [sol(1), sol(3), sol(5)];
B = [sol(2), sol(4), sol(6)];

for j = 1:N
    x = X(j);
    [Phi(j), index(j)] = phi_i(x,Omega,A,B,D,S);
end

I1 = find(index==1);
I2 = find(index==2);
I3 = find(index==3);

% Phin(I1) = (1/sqrt(3))*Phi(I1)/sqrt(trapz(X(I1),Phi(I1).^2));
% Phin(I2) = (1/sqrt(3))*Phi(I2)/sqrt(trapz(X(I2),Phi(I2).^2));
% Phin(I3) = (1/sqrt(3))*Phi(I3)/sqrt(trapz(X(I3),Phi(I3).^2));

Phin = Phi/sqrt(trapz(X,Phi.^2));
end

