function [elements2edges,faces,edges,boundaryFaces,boundaryNodes,boundaryEdges,boundaryEdgeIndex] = prepareMesh(elements)
    nE = size(elements,1);

    %Obtain Edges
    totalEdge = uint32([elements(:,[1 2]); elements(:,[1 3]); elements(:,[1 4]); ...
                        elements(:,[2 3]); elements(:,[2 4]); elements(:,[3 4])]);
    sortedTotalEdge = sort(totalEdge,2);
        [edges, ~, je] = unique(sortedTotalEdge,'rows');
    %Obtain Faces
    totalFace = uint32([elements(:,[2 3 4]); elements(:,[1 4 3]); ...
                        elements(:,[1 2 4]); elements(:,[1 3 2])]);
    sortedTotalFace = sort(totalFace,2);                
    [faces, ~, ~] = unique(sortedTotalFace,'rows');

    %Obtain Elements2edges
    elements2edges = uint32(reshape(je,nE,6));

    %Obtain Boundary Faces
    [~, boundaryFaces] = findboundary3(elements) ;

    % Obtain Boundary Edges 
    totalBoundaryEdges = uint32([boundaryFaces(:,[1,2]);boundaryFaces(:,[1,3]);boundaryFaces(:,[2,3])]);
    sortedBoundaryEdges = sort(totalBoundaryEdges,2);
    [boundaryEdges,~,~] = unique(sortedBoundaryEdges,'rows');
    
    % Obtain Boundary nodes
    boundaryNodes = unique(sort(boundaryFaces(:)));
    
    % Obtain Boundary edges Index
    [~, boundaryEdgeIndex]=ismember(boundaryEdges,edges,'rows');

end