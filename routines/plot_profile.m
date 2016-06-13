% Plots different observables at the final time. These include
%   - profile of the occupation number
%   - map of the correlation matrix
%   - profile of the correlation of the middle site with the others

figure
plot(1:N,real(n(:,end)))
xlabel('site $i$');
ylabel('$\langle n_i \rangle$');
ylim([0 Inf])
xlim([1 N])
% set(gca,'XTick',1:N)
title(sprintf('$N$=%.0d, $d$=%.0d, d$t$=%.0d, $\\epsilon$=%.0d',N,d,dt,epsilon));

figure
imagesc(real(g2(:,:,end)))
xlabel('site $i$');
ylabel('site $j$');
axis equal 

figure
plot(1:N,real(g2(ceil(N/2),:,end)))
xlim([1 N])
ylim([0 inf])
xlabel('site $j$');
ylabel('$g^{(2)}_{\diamond,j}$')