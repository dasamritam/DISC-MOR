function T = T_i0(x,X)

x_1 = X(1);
x_2 = X(2);
x_3 = X(3);
x_4 = X(4);

A = 10;
a = 10;

T_ambient = 20;

if (x_1 <= x) && (x <  x_2)
    T = 15;
    %T = a*sin(A*pi/6*x);
end

if (x_2 <= x) && (x <  x_3)
    T = 40;
    %T = a*sin(A*pi/6*x);
end

if (x_3 <= x) && (x <=  x_4)
    T = 15;
    %T = a*sin(A*pi/6*x);
end

T = T + T_ambient;

end

