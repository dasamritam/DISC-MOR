function G = g_i2(x,K,H,X)

K_1 = K(1);
K_2 = K(2);
K_3 = K(3);

h_1 = H(1);
h_2 = H(2);
h_3 = H(3);
h_4 = H(4);

x_1 = X(1);
x_2 = X(2);
x_3 = X(3);
x_4 = X(4);

if (x_1 <= x) && (x <  x_2)
    G = (K_2*K_3*h_2*h_3*h_4*(K_1 + h_1*x - h_1*x_1))/(K_1*K_2*K_3*h_1*h_2*h_3 + K_1*K_2*K_3*h_1*h_2*h_4 + K_1*K_2*K_3*h_1*h_3*h_4 + K_1*K_2*K_3*h_2*h_3*h_4 - K_1*K_2*h_1*h_2*h_3*h_4*x_3 - K_1*K_3*h_1*h_2*h_3*h_4*x_2 - K_2*K_3*h_1*h_2*h_3*h_4*x_1 + K_1*K_2*h_1*h_2*h_3*h_4*x_4 + K_1*K_3*h_1*h_2*h_3*h_4*x_3 + K_2*K_3*h_1*h_2*h_3*h_4*x_2);
 
end

if (x_2 <= x) && (x <  x_3)
    G = (K_3*h_3*h_4*(K_1*K_2*h_1 + K_1*K_2*h_2 + K_1*h_1*h_2*x - K_1*h_1*h_2*x_2 - K_2*h_1*h_2*x_1 + K_2*h_1*h_2*x_2))/(K_1*K_2*K_3*h_1*h_2*h_3 + K_1*K_2*K_3*h_1*h_2*h_4 + K_1*K_2*K_3*h_1*h_3*h_4 + K_1*K_2*K_3*h_2*h_3*h_4 - K_1*K_2*h_1*h_2*h_3*h_4*x_3 - K_1*K_3*h_1*h_2*h_3*h_4*x_2 - K_2*K_3*h_1*h_2*h_3*h_4*x_1 + K_1*K_2*h_1*h_2*h_3*h_4*x_4 + K_1*K_3*h_1*h_2*h_3*h_4*x_3 + K_2*K_3*h_1*h_2*h_3*h_4*x_2);
end

if (x_3 <= x) && (x <=  x_4)
   G = 1 - (K_1*K_2*K_3*h_1*h_2*h_3 - K_1*K_2*h_1*h_2*h_3*h_4*x + K_1*K_2*h_1*h_2*h_3*h_4*x_4)/(K_1*K_2*K_3*h_1*h_2*h_3 + K_1*K_2*K_3*h_1*h_2*h_4 + K_1*K_2*K_3*h_1*h_3*h_4 + K_1*K_2*K_3*h_2*h_3*h_4 - K_1*K_2*h_1*h_2*h_3*h_4*x_3 - K_1*K_3*h_1*h_2*h_3*h_4*x_2 - K_2*K_3*h_1*h_2*h_3*h_4*x_1 + K_1*K_2*h_1*h_2*h_3*h_4*x_4 + K_1*K_3*h_1*h_2*h_3*h_4*x_3 + K_2*K_3*h_1*h_2*h_3*h_4*x_2);
end
 

end
