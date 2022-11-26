%%
nc = size(Coordinates,1);
ne = size(Elements,1);
ned = size(edges,1);
% Get coordinates
X = Coordinates(:,1);
Y = Coordinates(:,2);
Z = Coordinates(:,3);
% Get Degrees of freedom 
u_h = x(1:dimU);


%% 
% Preprocessing data 
%CData_u = zeros(nc,1);
%for jdx=1:ne
%    CData_u(Elements(jdx,:)) = u_h(jdx);
%end 

%curlu_exact = curlu(coordinates);
u_exact = u(Coordinates); 
w_exact = w(Coordinates);
%%
% Plot 
figure
patch(...
    'Vertices',[X,Y,Z],...
    'Faces', faces,'CData',u_h,'FaceColor','interp','EdgeColor','black');
    %'Marker', plot_param.vertice_marker,...
    %'FaceColor', plot_param.face_color,...
    %'EdgeColor', plot_param.edge_color,...
    %'FaceAlpha', plot_param.face_alpha,...
    %'LineWidth', plot_param.edge_width...
    %);
view([45,45])
title('Solución exacta u1')

% transversal section

figure
patch(...
    'Vertices',[X,Y,Z],...
    'Faces', faces,'CData',u_h,'FaceColor','interp','EdgeColor','black');
    %'Marker', plot_param.vertice_marker,...
    %'FaceColor', plot_param.face_color,...
    %'EdgeColor', plot_param.edge_color,...
    %'FaceAlpha', plot_param.face_alpha,...
    %'LineWidth', plot_param.edge_width...
    %);
view([45,45])
title('Solución numerica u1')

