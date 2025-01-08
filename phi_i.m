function [phi, index] = phi_i(x,omega_,A,B,D,X)

D_1 = D(1);
D_2 = D(2);
D_3 = D(3);

x_1 = X(1);
x_2 = X(2);
x_3 = X(3);
x_4 = X(4);

if (x_1 <= x) && (x <  x_2)
    phi = A(1)*cos((omega_/D_1)*x) + B(1)*sin((omega_/D_1)*x);
    index = 1;
end

if (x_2 <= x) && (x <  x_3)
    phi = A(2)*cos((omega_/D_2)*x) + B(2)*sin((omega_/D_2)*x);
    index = 2;
end

if (x_3 <= x) && (x <=  x_4)
    phi = A(3)*cos((omega_/D_3)*x) + B(3)*sin((omega_/D_3)*x);
    index = 3;
end

end

