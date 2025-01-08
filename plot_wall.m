function  plot_wall(X,Tmax,Tmin)

cla;

color{1} = [0.6 0.6 0.6];
color{2} = [0.8 0.8 0.8];
color{3} = [0.7 0.7 0.7];

for i = 1:length(X)-1
    lay{i}.vertices = [X(i) Tmin; X(i+1) Tmin; X(i) Tmax; X(i+1) Tmax];
    lay{i}.faces = [1 2 4 3];
end

for i = 1:length(X)-1
    patch(lay{i},'facealpha',0.8,'facecolor',color{i},'linestyle','none');
end

hold on;

for i = 2:length(X)-1
    line = [X(i)*ones(5,1) , linspace(Tmin,Tmax,5).'];
    plot(line(:,1),line(:,2),'k-.');
end

%xlabel('Postition $x$ (cm)','interpreter','latex');
%ylabel('Temperature $T(x,t)$ ($^\circ$C)','interpreter','latex');
grid on; box off;

end

