function [A,I,J] = buildAMatLap(coordinates,elements,elem2edges,edges,faces,...
    pU)

nE = size(elements,1);
nF = size(faces,1);
nEd = size(edges,1);
nC = size(coordinates,1);


A = [];
I = []; 
J = [];

dimU = nC+(pU-1)*nEd; 

%% (grad u,grad w)
[Atmp,Itmp,Jtmp] = buildBLFGradUGradW(coordinates,elements,elem2edges,edges,pU);

I = [I;Itmp];
J = [J;Jtmp];
A = [A;Atmp];

%% OUTPUT
if(nargout==1)
    A = sparse(I,J,A,dimU,dimU);
end