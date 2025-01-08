function [A,symbols] = construct_GAMMA()

syms a1 b1 a2 b2 a3 b3 a4 b4

syms x x_1 x_2 x_3 x_4 h_1 h_2 h_3 h_4 K_1 K_2 K_3 K_4 D_1 D_2 D_3 D_4 omega

M = 3;

a{1} = a1;
a{2} = a2;
a{3} = a3;

b{1} = b1;
b{2} = b2;
b{3} = b3;

K{1} = K_1; K{2} = K_2; K{3} = K_3; K{4} = K_4;
h{1} = h_1; h{2} = h_2; h{3} = h_3; h{4} = h_4;
D{1} = D_1; D{2} = D_2; D{3} = D_3; D{4} = D_4;
X{1} = x_1; X{2} = x_2; X{3} = x_3; X{4} = x_4;

%% form general solution for g_i,1(x), and g_i,x(x)

for i = 1:M
    U{i} = a{i}*cos((omega/D{i})*x) + b{i}*sin((omega/D{i})*x);
end

%% set boundaries for U_i(x_i)

% surface boundary conditions
f{1} = K{1}*subs(diff(U{1},x),x,X{1}) - h{1}*subs(U{1},x,X{1}); 
f{2} = K{M}*subs(diff(U{M},x),x,X{M+1}) + h{M+1}*subs(U{M},x,X{M+1});

% internal boundary conditions
for i = 1:M-1
    f{i+2} = -K{i}*subs(diff(U{i},x),x,X{i+1}) - h{i+1}*(subs(U{i},x,X{i+1}) - subs(U{i+1},x,X{i+1}));
end

for i = 1:M-1
    f{i+4} = K{i}*subs(diff(U{i},x),x,X{i+1}) - K{i+1}*subs(diff(U{i+1},x),x,X{i+1});
end

% form the vector of unknowns
c = [a{1},b{1},a{2},b{2},a{3},b{3}].';

for i = 1:6
   A(:,i) = jacobian(f{i},c);
end

A = A.';

symbols = [x_1 x_2 x_3 x_4 h_1 h_2 h_3 h_4 K_1 K_2 K_3 D_1 D_2 D_3 omega];
end

