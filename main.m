

figure

I_s = [0, 20, 34,  45];

for ii = 1:4
sol = neuron_model(600, 0.5, I_s(ii), 0, 50, 550);
subplot(4,2, 2*(ii-1)+1)
plot(sol.t, sol.y(:, 1));
ylabel('V_s')
axis([0 600 -70 70])
title(['Current to soma is: ', num2str(I_s(ii)), 'uA']);

subplot(4,2, 2*(ii-1)+2)
plot(sol.t, sol.y(:, 3));
ylabel('V_d')
axis([0 600 -70 0])


end

