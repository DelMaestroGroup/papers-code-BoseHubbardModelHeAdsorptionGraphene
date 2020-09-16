clear all

Wanner_3D_x = readmatrix("../../Data/Wanner_3D_x.txt");
Wanner_3D_y = readmatrix("../../Data/Wanner_3D_y.txt");
Wanner_3D_rho = readmatrix("../../Data/Wanner_3D_rho.txt");

ind_shift = 7;
hex_x_ind1 = 26-ind_shift:48:577;
hex_y_ind1 = 26-ind_shift:48:577;

for i = 1:12
   for j = 1:12
       hex_x1(i, j) = Wanner_3D_x(hex_x_ind1(i), hex_y_ind1(j));
       hex_y1(i, j) = Wanner_3D_y(hex_x_ind1(i), hex_y_ind1(j));
       hex_z1(i, j) = 0;%Wanner_3D_rho(hex_x_ind1(i), hex_y_ind1(j))+0.02;
   end
end

hex_x_arr1 = reshape(hex_x1, 1, 144);
hex_y_arr1 = reshape(hex_y1, 1, 144);
hex_z_arr1 = reshape(hex_z1, 1, 144);

hex_x_ind2 = 26+ind_shift:48:577;
hex_y_ind2 = 26+ind_shift:48:577;

for i = 1:12
   for j = 1:12
       hex_x2(i, j) = Wanner_3D_x(hex_x_ind2(i), hex_y_ind2(j));
       hex_y2(i, j) = Wanner_3D_y(hex_x_ind2(i), hex_y_ind2(j));
       hex_z2(i, j) = 0;%Wanner_3D_rho(hex_x_ind2(i), hex_y_ind2(j))+0.02;
   end
end

hex_x_arr2 = reshape(hex_x2, 1, 144);
hex_y_arr2 = reshape(hex_y2, 1, 144);
hex_z_arr2 = reshape(hex_z2, 1, 144);

figure(1)
mesh(Wanner_3D_x, Wanner_3D_y, Wanner_3D_rho, 'EdgeAlpha', 0.4, 'FaceAlpha', 0)
hold on
scatter3(hex_x_arr1, hex_y_arr1, hex_z_arr1, 'filled', 'MarkerFaceColor', 'black');
hold on
scatter3(hex_x_arr2, hex_y_arr2, hex_z_arr2, 'filled', 'MarkerFaceColor', 'black');
hold off
set(gca, 'FontName', 'Computer Modern Roman', 'FontSize', 14)
xlim([-8,8])
ylim([-8,8])
xlabel('x(\AA)', 'interpreter','latex','FontName','Computer Modern Roman', ...
       'fontsize', 16)
ylabel('y(\AA)', 'interpreter','latex','FontName','Computer Modern Roman', ...
       'fontsize', 16)
zlabel('$\rho$', 'interpreter','latex','FontName','Computer Modern Roman', ...
       'fontsize', 16, 'Rotation', 0)
xticks([-8, -4, 0, 4, 8])
xticklabels({'-8','-4','0', '4','8'})
yticks([-8, -4, 0, 4, 8])
yticklabels({'-8','-4','0', '4','8'})

saveas(figure(1), "Wannier3D.svg")