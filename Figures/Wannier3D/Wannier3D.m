clear all

Wanner_3D_x = readmatrix("../../Data/Wanner_3D_x.txt");
Wanner_3D_y = readmatrix("../../Data/Wanner_3D_y.txt");
Wanner_3D_rho = readmatrix("../../Data/Wanner_3D_rho.txt");

ind_shift = 8;
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
hold on 

for i1 = 1:12
    for j1 = 1:12
        for i2 = i1-2:i1+2
            for j2 = j1-2:j1+2
                if i2 > 0 && i2 < 13 && j2 > 0 && j2 < 13
                    if sqrt((hex_x1(i1, j1) - hex_x2(i2, j2))^2 + ...
                            (hex_y1(i1, j1) - hex_y2(i2, j2))^2) < 2
                        plot3([hex_x1(i1, j1), hex_x2(i2, j2)], ...
                              [hex_y1(i1, j1), hex_y2(i2, j2)], ...
                              [0, 0], 'Color', 'black', 'LineWidth', 1.5)
                        hold on
                    end
                end
            end
        end
    end
end

hold off

set(gca, 'FontName', 'Computer Modern Roman', 'FontSize', 14)
xlim([-6,6])
ylim([-6,6])
xlabel('x(\AA)', 'interpreter','latex','FontName','Computer Modern Roman', ...
       'fontsize', 16)
ylabel('y(\AA)', 'interpreter','latex','FontName','Computer Modern Roman', ...
       'fontsize', 16)
zlabel('$\rho$', 'interpreter','latex','FontName','Computer Modern Roman', ...
       'fontsize', 16, 'Rotation', 0)
xticks([-6, -3, 0, 3, 6])
xticklabels({'-6','-3','0', '3','6'})
yticks([-6, -3, 0, 3, 6])
yticklabels({'-6','-3','0', '3','6'})
%axis('square')

saveas(figure(1), "Wannier3D.svg")