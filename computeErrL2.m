function err = computeErrL2(coordinates,elements,elements2edges,edges,uh,u,pU)
    %dim = size(uh,1);
    %dimP = 18;
    nE = size(elements,1);
    nC = size(coordinates,1);
    nEd = size(edges,1);
    err = zeros(nE,1);
    %dimP = 4*pU+6*(pU-1); 
    %dimU = nC+(pU-1)*nEd;
    %*** Quadrature rule
    [wt,xa,ya,za]=QuadratureRule(1);
    %ia = [1,2,1,1];
    %nQ = length(wt);
    %*** Volumes 
    d21 = coordinates(elements(:,2),:)-coordinates(elements(:,1),:);
    d31 = coordinates(elements(:,3),:)-coordinates(elements(:,1),:);
    d41 = coordinates(elements(:,4),:)-coordinates(elements(:,1),:);
    vol2 = d21(:,1).*(d31(:,2).*d41(:,3)-d31(:,3).*d41(:,2))...
    -d31(:,1).*(d21(:,2).*d41(:,3)-d21(:,3).*d41(:,2))...
    +d41(:,1).*(d21(:,2).*d31(:,3)-d21(:,3).*d31(:,2));
    vol2 = abs(vol2);
    
    for jdx=1:nE
        % ***Forma 1 
        p = coordinates(elements(jdx,:),:);
        F = [d21(jdx,:)',d31(jdx,:)',d41(jdx,:)']*[xa;ya;za]+coordinates(elements(jdx,1),:)';
        U = u([F(1,:)',F(2,:)',F(3,:)'])';
        [phi,~,~,~] = basis_u(F(1,:),F(2,:),F(3,:),p');
        if(pU > 0)
            Itmp = elements(jdx,:);
        end
        if(pU>1)
            Itmp = [elements(jdx,:),elements2edges(jdx,:)+nC];  
        end
        errorT = (uh(Itmp)'*phi-U(1,:)).^2*wt';
        err(jdx,1) = errorT*vol2(jdx);  
    end
end