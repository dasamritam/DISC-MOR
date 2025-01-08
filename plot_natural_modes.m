function plot_natural_modes(PHI,OMEGA,X,S,winx,winy,style)

N = winx*winy;

for i = 1:N
subplot(winx,winy,i);
plot_wall(S,2,-2); 
hold on;
plot(X,PHI(i,:),style);
text = ['$\omega_{', num2str(i),'}=$',num2str(OMEGA(i),'%10.2f\n'),' (rad/s)'];
    title(text,'interpreter','latex','fontsize',12);
axis equal
axis tight
end

