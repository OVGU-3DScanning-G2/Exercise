% data =load('cone.xyz');
% data =load('sphere.xyz');
data =load('sphere2.xyz');
data_size = size(data)




close all
hold on;
axis equal;
% title ('sphere of radius ');
X = data(:,1);
Y = data(:,2);
Z = data(:,3);
scatter3(X,Y,Z)
disp('sphere_fit by Sumith YD, Syracuse University')
tic
[Center,Radius] = sphere_fit(X,Y,Z)
toc
[x, y, z] = sphere (40);
 
surf (Radius*x+Center(1), Radius*y+Center(2), Radius*z+Center(3))%,'EdgeColor', 'none','interp','FaceAlph‌​a',0.1)
% set(s,'facealpha',0.5);
% set(s,'edgecolor','none'); 

% alpha (0.99);
disp('Function by Alan Jennings, University of Dayton')
tic
[Center1,Radius1] = sphereFit(data);
Center1'
Radius1
toc
[x1, y1, z1] = sphere (40);
% alpha 0.5
% surf (Radius1*x1+Center1(1), Radius1*y1+Center1(2), Radius1*z1+Center1(3))
euclidean_difference_of_centers = norm(Center - Center1',2)
difference_in_radius = Radius - Radius1