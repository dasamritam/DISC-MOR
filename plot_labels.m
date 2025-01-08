function plot_labels(ypos, max_temp, S)

delta = 0.125;

axis([min(S) max(S) 0 max_temp])
text((S(2)-S(1))/2 + S(1)-delta,ypos,'Layer 1','interpreter','latex','fontsize',12)
text((S(3)-S(2))/2 + S(2)-delta,ypos,'Layer 2','interpreter','latex','fontsize',12)
text((S(4)-S(3))/2 + S(3)-delta,ypos,'Layer 3','interpreter','latex','fontsize',12)

end

