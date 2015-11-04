%rigid_transform_2D  Returns optimal rigid/Euclidean transform in 2D
%
%   Returns R, a 2x2 rotation matrix and T, a
%   2x1 translation vector which transforms the points in B to A using
%   single value decomposition. Point matrices A and B are Nx3, with 3x3
%   being the smallest size.
%
%   Example: [R,T]=rigid_transform(A,B) will give R and T so
%   Anew=R*A'+repmat(T,1,size(A,1)) will provide Anew' which is
%   approximately equal to B.
%                    
%
%   See also SVD.
%   
%   Copyright 2015 M. J. Roy
%   $Revision: 1.0$  $Date: 2015/10/30$
function [R,t] = rigid_transform_2D(A, B)
    if nargin ~= 2
	    fprintf('Missing parameters in rigid_transform\n');
    end

    assert(isequal(size(A),size(B)),...
        'Number of points to compare isn''t the same');
    assert(size(A,1)>=3,...
        'Must have more than 3 points to compare.');
    
    centroid_A = mean(A);
    centroid_B = mean(B);

    N = size(A,1);

    H = (A - repmat(centroid_A, N, 1))' * (B - repmat(centroid_B, N, 1));

    [U,~,V] = svd(H);

    R = V*U';

    if det(R) < 0
        fprintf('Reflection detected in rigid_transform\n');
        V(:,2) = -1;
        R = V*U';
    end

    t = -R*centroid_A' + centroid_B';
end