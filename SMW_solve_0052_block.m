function  X = SMW_solve_0052_block(A,b,M,P,Q,blk)    
    [~,n] = size(P);
    v = Q;     
    m1 = blk;
    p = n/m1;
    Is = eye(m1); 
    X = M \ b;
    Y = M \ P;

    for k = 1:p-1
        I = ((k-1)*m1+1):k*m1';       
        sI=((k-1)*m1+1);
        eI=k*m1;
        Xblock=X(sI:eI);
        Yblock=Y(sI:eI,sI:eI);
        %%% previous variant : X = X - Y(:,I) * (((Is + v(:,I)' * Y(:,I))) \ (v(:,I)' * X));
        % v(:,I)' * Y(:,I) and (v(:,I)' * X) - v is identity matrix.This matrix multiplication
        % is not necessary, because, in essence, this is the same thing as choosing a sub-matrix with
        % the row and columns with indices I.
        % i.e Xblock=X(I)=X(sI:eI)=v(:,I)' * X  and Yblock=Y(I,I)=Y(sI:eI,sI:eI);= v(:,I)' * Y(:,I)
        X = X - Y(:,I) * vector_update(Xblock,Yblock); 
        for j = k+1:p
            J = ((j-1)*m1+1):j*m1;
            sJ=((j-1)*m1+1);
            eJ=j*m1;
            Xblock=Y(sI:eI,sJ:eJ);
            Yblock=Y(sI:eI,sI:eI);        
            Y(:,J) = Y(:,J) - Y(:,I) * vector_update(Xblock,Yblock); 
        end
    end
 
    k = k+1;
    I = ((k-1)*m1+1):k*m1;    
    sI=((k-1)*m1+1);
    eI=k*m1;     
    Xblock=X(sI:eI);
    Yblock=Y(sI:eI,sI:eI);
    X = X - Y(:,I) * vector_update(Xblock,Yblock);        
end

function X=vector_update(Xblock,Yblock)    
    Is=eye(size(Yblock));
    X=(Is+Yblock) \ Xblock;    
end