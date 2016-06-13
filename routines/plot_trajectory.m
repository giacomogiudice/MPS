% Plots the trajectory of the computed observables after a simulation as a
% function of time. These include (if available): 
%   - occupation number for each site
%   - same-site correlation
%   - nearest neighbor correlation
%   - largest bond dimension (for MPS only)

figure
hold on
for i = 1:N
    plot(t_sampling,real(n(i,:)));
end
xlabel('$t$');
ylabel('$\langle n \rangle$');
title(sprintf('$N$=%.0d, $d$=%.0d, d$t$=%.0d',N,d,dt));

if exist('g2','var')
    figure
    hold on
    for i = 1:N
        plot(t_sampling,squeeze(real(g2(i,i,:))));
    end
    xlabel('$t$');
    ylabel('$g^{(2)}_{i,i}$');
    title(sprintf('$N$=%.0d, $d$=%.0d, d$t$=%.0d',N,d,dt));
    xlim([0 T])
    figure
    hold on
    for i = 1:N-1
        plot(t_sampling,squeeze(real(g2(i,i+1,:))));
    end
    xlabel('$t$');
    ylabel('$g^{(2)}_{i,i+1}$');
    title(sprintf('$N$=%.0d, $d$=%.0d, d$t$=%.0d',N,d,dt));
    xlim([0 T])   
end

if exist('bond_size','var')
    figure
    plot(t_sampling,bond_size,'.')
    xlabel('$t$');
    ylabel('$\max{D_b}$');
    ylim([0 D_max]) 
end