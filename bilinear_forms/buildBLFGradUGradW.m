function [A,I,J] = buildBLFGradUGradW(coordinates,elements,elem2edges,edges,pU)
nC = size(coordinates,1);
nE = size(elements,1);
nEd = size(edges,1);
dimP = 4*pU+6*(pU-1); 
% building matrix
A = zeros(dimP*dimP*nE,1);
I = zeros(dimP*dimP*nE,1);
J = zeros(dimP*dimP*nE,1);

% volumes 
d21 = coordinates(elements(:,2),:)-coordinates(elements(:,1),:);
d31 = coordinates(elements(:,3),:)-coordinates(elements(:,1),:);
d41 = coordinates(elements(:,4),:)-coordinates(elements(:,1),:);
vol2 = d21(:,1).*(d31(:,2).*d41(:,3)-d31(:,3).*d41(:,2))...
    -d31(:,1).*(d21(:,2).*d41(:,3)-d21(:,3).*d41(:,2))...
    +d41(:,1).*(d21(:,2).*d31(:,3)-d21(:,3).*d31(:,2));
vol2 = abs(vol2);
% quadrature rule 
[wt,xa,ya,za] = QuadratureRule(1);
for jdx=1:nE
    p = coordinates(elements(jdx,:),:);
    F = [d21(jdx,:)',d31(jdx,:)',d41(jdx,:)']*[xa;ya;za]+coordinates(elements(jdx,1),:)';
    [~,dxu,dyu,dzu]=basis_u(F(1,:),F(2,:),F(3,:),p',pU);
    % Calculate a(u,v) over K1,K2,K3 and K4
    gradugradv = repelem(dxu,dimP,1).*repmat(dxu,dimP,1)...
            +repelem(dyu,dimP,1).*repmat(dyu,dimP,1)...
            +repelem(dzu,dimP,1).*repmat(dzu,dimP,1);
           
    gradugradv = gradugradv*wt';
    A((jdx-1)*dimP*dimP+1:jdx*dimP*dimP) = gradugradv';
    if(pU > 0)
        Itmp = repelem(elements(jdx,:),dimP,1)';
        Jtmp = repmat(elements(jdx,:),dimP,1);
    end
    if(pU>1)
        Itmp = repelem([elements(jdx,:),elem2edges(jdx,:)+nC],dimP,1)';  
        Jtmp = repmat([elements(jdx,:),elem2edges(jdx,:)+nC],dimP,1);
    end
    I((jdx-1)*dimP*dimP+1:jdx*dimP*dimP) = Itmp(:);
    J((jdx-1)*dimP*dimP+1:jdx*dimP*dimP) = Jtmp(:);
end
end