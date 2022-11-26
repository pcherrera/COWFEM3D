%% %% Poisson Problem FEM 3D
% Weak formulation 
% 
% MODEL PROBLEM
% % - div grad w = f    in Omega
% %            u = o   on Gamma 

%% 
clear all
close all
addpath bilinear_forms/
%% Create a mesh 
%[Coordinates, Elements] = cubemesh([0,1,0,1,0,1]);
load('Vaca.mat');
% obtain boundary faces
%h = 1;
%coordinates = [0,0,0;h,0,0;0,h,0;0,0,h];
%elements = [1,2,3,4];

%% Set Polynomial orders 
pU = 1;
markedE = [];
maxL = 4;
nEmax = 4000;

%% Set adaptivity parameter \theta
% 0<theta<1 is adaptive refinement using 
% theta = 1 is uniform mesh refinement
% theta = 1;
%% example 1 smooth function in (0,1).^3
u = @(x) x(:,1).*x(:,2).*x(:,3).*(x(:,1)-1).*(x(:,2)-1).*(x(:,3)-1);
w = @(x) [x(:,1).*x(:,2).*x(:,3).*(x(:,2)-1).*(x(:,3)-1)+x(:,2).*x(:,3).*(x(:,1)-1).*(x(:,2)-1).*(x(:,3)-1),...
        x(:,1).*x(:,2).*x(:,3).*(x(:,1)-1).*(x(:,3)-1)+x(:,1).*x(:,3).*(x(:,1)-1).*(x(:,2)-1).*(x(:,3)-1),...
        x(:,1).*x(:,2).*x(:,3).*(x(:,1)-1).*(x(:,2)-1)+x(:,1).*x(:,2).*(x(:,1)-1).*(x(:,2)-1).*(x(:,3)-1)];
f = @(x) -2.*x(:,1).*x(:,2).*(x(:,1) - 1).*(x(:,2) - 1) - 2*x(:,1).*x(:,3).*(x(:,1) - 1).*(x(:,3) - 1) - 2*x(:,2).*x(:,3).*(x(:,2) - 1).*(x(:,3) - 1);
g = @(x) u(x);
for j=1:1
    %refine
    %[Coordinates,Elements] = uniformrefine3(Coordinates,Elements);
    
    %Obtain geometric data structures 
    [elements2edges,faces,edges...
        ,boundaryFaces,boundaryNodes,boundaryEdges,boundaryEdgeIndex] = prepareMesh(Elements);
    % store some values(numbers of elements, vertices,etc.)
    nE(j) = size(Elements,1);
    nC(j) = size(Coordinates,1);
    nEd(j) = size(edges,1);
    fprintf('nE: %d, nC: %d, nEd: %d. \n',nE(j),nC(j),nEd(j));
    % dimension of spaces
    dimU = nC(j)+(pU-1)*nEd(j);
    % degrees of freedom 
    DOF(j) = dimU;
    % build bilinear form a(u,v)
    A = buildAMatLap(Coordinates,Elements,elements2edges,edges,faces,pU);
    %*** build RHS
    bVec = buildRHSLap(Coordinates,Elements,elements2edges,edges,f,pU);
    %build linear system
    dof = size(A,1);
    %jdx = [];
    jdx=boundaryNodes;
    %*** Boundary Nodes dof (Just for a cube domain)
    x = zeros(dof,1);
    x(boundaryNodes) = g(Coordinates(boundaryNodes,:));
    freenodes = setdiff(1:dof,jdx);
    ADir = A(freenodes,boundaryNodes)*x(boundaryNodes);
    %*** Solve A*xfree = F-Adir*xDir
    x(freenodes) = A(freenodes,freenodes)\(bVec(freenodes)-ADir);
    %*** Residual error estimator
    %*** L2 error of u
    tmpErr = computeErrL2(Coordinates,Elements,elements2edges,edges,x,u,pU);
    errU(j) = sqrt(sum(tmpErr));
    %*** L2 error of curlcurlu
    %tmpErr1 = computeErrL2_c(coordinates,elements,elements2edges,x,w);
    %errW(j) = sqrt(sum(tmpErr1));
    %err(j) = sqrt(x'*bVec-0.5*x'*B*x);
    %unorm(j)=x'*A*x;
end

%% Plot convergence rates
%close all
minP=min(pU)+1;
%minP = 2/3;
figure
loglog(DOF,errU,'--o',...
    DOF,errU(1)*(DOF/DOF(1)).^(-(minP)/3),'--k','LineWidth',2.5);
legend('err eh');
xlabel('grados de libertad')
ylabel('Error L2')