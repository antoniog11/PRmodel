% replicate to figure 1b
I_s = [20, 34,  45];
I_s = I_s(end:-1:1);

figure
set(gcf, 'position', [1 1 1200 800])
for ii = 1:3
sol = neuron_model(600, 3, I_s(ii), 0, 50, 550);
subplot(3,1, ii)
plot(sol.t, sol.y(:, 1), 'LineWidth', 1.5);
ylabel('V_s(mV)')
axis([0 600 -80 50])
title(['I_s=', num2str(I_s(ii)), 'uA/cm^2']);
box off
ax = gca;
ax.FontSize = 15;
ax.LineWidth = 2;

end
xlabel('t(ms)')


%%
g_c = 0.9;
sol = neuron_model(600, g_c, 0, 75, 1, 600);

figure
set(gcf, "Position", [1 1 500 1000])
subplot 411
plot(sol.t, sol.y(:, 1))
box off
ylabel('V_S(mV)')
title('g_{Ca} = 0 mS/cm^2')
axis([0 600 -80 50])
ax=gca;
ax.FontSize=15;
ax.LineWidth = 1.5;

subplot 412
plot(sol.t, sol.y(:, 3))
box off
ylabel('V_D(mV)')
% axis([0 600 -80 120])
axis([0 600 -80 20])
ax=gca;
ax.FontSize=15;
ax.LineWidth = 1.5;


subplot 413
plot(sol.t, sol.y(:,4).*sol.y(:, 5).*0.*(sol.y(:, 3) - 120));
box off
ylabel('I_{Ca}(uA/cm^2)')
% axis([0 600 -550 10])
axis([0 600 -80 50])
ax=gca;
ax.FontSize=15;
ax.LineWidth = 1.5;

subplot 414
plot(sol.t, g_c*(sol.y(:, 1)-sol.y(:,3)));
box off
ylabel('I_{DS}(uA/cm^2)')
% axis([0 600 -160 50])
axis([0 600 -80 50])
ax=gca;
ax.FontSize=15;
ax.LineWidth = 1.5;

%% test
g_c =1e-1;
sol = neuron_model(600, g_c, 0, 100, 0, 600);

figure
subplot 311
plot(sol.t, sol.y(:, 1))
subplot 312
plot(sol.t, sol.y(:, 3))
subplot 313
plot(sol.t, g_c*(sol.y(:, 1)-sol.y(:,3)));

figure
subplot 311
plot(sol.t, sol.y(:, 4))
subplot 312
plot(sol.t, sol.y(:, 5))
subplot 313
plot(sol.t, 40*sol.y(:, 5).*sol.y(:, 4))

figure
plot(sol.t, sol.y(:,4).*sol.y(:, 5).*40.*(sol.y(:, 3) - 120));