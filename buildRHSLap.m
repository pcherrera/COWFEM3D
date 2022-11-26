% NO TOCAR 
function bVec = buildRHSLap(coordinates,elements,elements2edges,edges,f,pU)
    nC = size(coordinates,1);
    nE = size(elements,1);
    nEd = size(edges,1);
    dimP = 4*pU+6*(pU-1); 
    dimU = nC+(pU-1)*nEd;
    % building matrix
    A = zeros(dimP*nE,1);
    I = zeros(dimP*nE,1);
    J = ones(dimP*nE,1);

    % Quadrature points
    [wt,xa,ya,za]=QuadratureRule(1);
    % Volumes 
    d21 = coordinates(elements(:,2),:)-coordinates(elements(:,1),:);
    d31 = coordinates(elements(:,3),:)-coordinates(elements(:,1),:);
    d41 = coordinates(elements(:,4),:)-coordinates(elements(:,1),:);
    vol2 = d21(:,1).*(d31(:,2).*d41(:,3)-d31(:,3).*d41(:,2))...
        -d31(:,1).*(d21(:,2).*d41(:,3)-d21(:,3).*d41(:,2))...
        +d41(:,1).*(d21(:,2).*d31(:,3)-d21(:,3).*d31(:,2));
    vol2 = abs(vol2);
    for jdx=1:nE
      p = coordinates(elements(jdx,:),:);
      F = [d21(jdx,:)',d31(jdx,:)',d41(jdx,:)']*[xa;ya;za]+coordinates(elements(jdx,1),:)';
      fT =  feval(f,[F(1,:)',F(2,:)',F(3,:)']);
      [u,~,~,~]=basis_u(F(1,:),F(2,:),F(3,:),p',pU);
      uf = (u*(fT(:,1).*wt')); 
      A((jdx-1)*dimP+1:jdx*dimP) = uf;
      if(pU > 0)
        Itmp = elements(jdx,:)';
      end
      if(pU > 1)
        Itmp = [elements(jdx,:),elements2edges(jdx,:)+nC]';  
      end
        I((jdx-1)*dimP+1:jdx*dimP) = Itmp(:);
    end 
    %% OUTPUT
    if(nargout==1)
        bVec = sparse(I,J,A,dimU,1);
    end
end