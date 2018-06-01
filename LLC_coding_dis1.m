function [Coeff] = LLC_coding_dis1(B, X,pos_mat,PB,knn, beta)

% ========================================================================
% USAGE: [Coeff]=LLC_coding_appr(B,X,knn,lambda)
% Approximated Locality-constraint Linear Coding
%
% Inputs
%       B       -M x d codebook, M entries in a d-dim space
%       X       -N x d matrix, N data points in a d-dim space
%       knn     -number of nearest neighboring
%       lambda  -regulerization to improve condition
%
% Outputs
%       Coeff   -N x M matrix, each row is a code for corresponding X
%
% Jinjun Wang, march 19, 2010
% ========================================================================
% function [Coeff] = LLC_coding_appr(B, X, knn, beta)

if ~exist('knn', 'var') || isempty(knn),
    knn = 5;
end

if ~exist('beta', 'var') || isempty(beta),
    beta = 1e-4;
end

nframe=size(X,1);
nbase=size(B,1);

% find k nearest neighbors
PX=pos_mat;
%PB=pos_mat(Backsuper_seq,:);
D=zeros(nframe,nbase);
for i=1:nframe
    for j=1:nbase
       D(i,j)=(PX(i,1)-PB(j,1)).^2+(PX(i,2)-PB(j,2)).^2;
    end
end
IDX = zeros(nframe, knn);
for i = 1:nframe,
	d = D(i,:);
	[dummy, idx] = sort(d, 'ascend');
	IDX(i, :) = idx(1:knn);
end

% llc approximation coding
II = eye(knn, knn);
Coeff = zeros(nframe, nbase);
for i=1:nframe
   idx = IDX(i,:);
   z = B(idx,:) - repmat(X(i,:), knn, 1);           % shift ith pt to origin
   C = z*z';                                        % local covariance
   C = C + II*beta*trace(C);                        % regularlization (K>D)
   w = C\ones(knn,1);
   w = w/sum(w);                                    % enforce sum(w)=1
   Coeff(i,idx) = w';
end
end
