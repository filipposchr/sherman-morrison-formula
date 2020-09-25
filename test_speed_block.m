%Test gia veltistopiimeni ekdosi tis SMW_solve_0052_block
for (i=1:100) 
    A=rand(2000,2000);
    n = 2000;
    b = 1:2000;
    b = b';
    M = diag(diag(A));
    P = A - M;
    Q = eye(n);
    blk = 20;
    tic; 
    [X]=SMW_solve_0052_block1(A,b,M,P,Q,blk);
    t1(i)=toc;
    tic;
    [X1]=SMW_solve_0052_block(A,b,M,P,Q,blk); 
    t2(i)=toc;
    nX(i)=norm(X-X1);
end
norm(nX)
sum(t1/t2)

%KWDIKAS SMW_solve_0052_block1 (mi veltistopiimeni sinartisi)
% function X = SMW_solve_0052_block1(A,b,M,P,Q,blk,sdir)
%     %block dimensions of n x blk
%     [~,n] = size(A);
%     m1 = blk; 
%     p = n/m1; %number of blocks
%     Is = eye(m1); %identity matrix of size blk
%     X = M \ b;
%     Y = M \ P;
% 
%     for k = 1:p-1 %go through all the blocks except the last one
%         I = ((k-1)*m1+1):k*m1'; %calculating indexing of blocks
%         X = X - Y(:,I) * (((Is + Q(:,I)' * Y(:,I))) \ (Q(:,I)' * X));
%      
%         for j = k+1:p %go through all the blocks
%             J = ((j-1)*m1+1):j*m1;
%             Y(:,J) = Y(:,J) - Y(:,I) * ((Is + Q(:,I)' * Y(:,I)) \ (Q(:,I)' * Y(:,J)));
%         end
%     end
%     %for the last block, compute X
%     k = k+1;
%     I = ((k-1)*m1+1):k*m1;
%     X = X - Y(:,I) * (((Is + Q(:,I)' * Y(:,I))) \ (Q(:,I)' * X));
% end