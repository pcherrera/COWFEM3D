%% GENERAR MALLA VAQUITA
model = createpde;
importGeometry(model,"Vaca.stl");
pdegplot(model);
mesh = generateMesh(model,'GeometricOrder','linear');
figure
pdemesh(model)
Coordinates = double(mesh.Nodes');
Elements = mesh.Elements';
save('Vaca','Coordinates','Elements');