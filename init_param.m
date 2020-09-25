function [A, b, M, P, Q, x_sol] = init_param(n,mx_id,sdir)
    switch sdir
        case 'colwise'
            A = MxMake_0052(mx_id,n);
            if (strcmp(mx_id,'wathen'))
                n = 3*n*11 + 2*n + 2*11 +1; %N = 3*NX*NY + 2*NX + 2*NY + 1
            end
            M = diag(diag(A));
            P = A - M;
            Q = eye(n);
        case 'rowwise' 
            A = MxMake_0052(mx_id,n);
            if (strcmp(mx_id,'wathen'))
                n = 3*n*11 + 2*n + 2*11 +1; %N = 3*NX*NY + 2*NX + 2*NY + 1
            end
            M = diag(diag(A));
            P = eye(n);
            Q = (A-M)';
        otherwise
            M = input('Enter an array M, values in []:')
            P = input('Enter an array P, values in []:')
            Q = eye(n);
            A = M + P*Q'
    end     
    m = size(A);
    n = m(1);
    c=1;
    x_sol = zeros(n,1);
    for i = 1:n
        if (mod(i,2) ~= 0)
            x_sol(i) = 1;
        else
            x_sol(i) = (-1)^(c+1) * (1/(2*i));
            c=c+1;
        end 
    end
    b = A * x_sol;
    
    cond_A = condest(A,1);
    testX = A \ b;
    X = SMW_solve_0052(A,b,M,P,Q,sdir)
    
    if (sdir == 'colwise')   
        exact_forward_error_colwise = zeros(n,1);
        sinepagomeno_forward_error_colwise = zeros(n,1);
        posteriori_backward_error_colwise = zeros(n,1);
  
        posteriori_backward_error_colwise =  norm(b-A*X,1) / (norm(A,1) * norm(X,1) * norm(b,1));
        sinepagomeno_forward_error_colwise = (2*cond_A) * norm(b-A*X,1) / (norm(A,1) * norm(X,1) + norm(b,1));
        exact_forward_error_colwise = norm(X - x_sol,1) / norm(x_sol,1);
            
        fprintf('----Sfalmata gia tin SMW_solve_0052,colwise, gia to mitroo %s-----', mx_id);
        fprintf('\nFragma Dikti Katastasis: %d\n',cond_A);
        fprintf('Sinepagwmeno Empros Sfalma: %d\n',sinepagomeno_forward_error_colwise);
        fprintf('A Posteriori Sxetiko Pisw Sfalma: %d\n',posteriori_backward_error_colwise);
        fprintf('Akrives sxetiko Pisw Sfalma: %d\n', exact_forward_error_colwise);
        fprintf('\n');
    end
    
    exact_forward_error_testX = zeros(n,1); 
    sinepagomeno_forward_error_testX = zeros(n,1); 
    posteriori_backward_error_testX = zeros(n,1);
    sinepagomeno_forward_error_testX = (2*cond_A) * norm(b-A*testX,1) / (norm(A,1) * norm(testX,1)+norm(b,1));
    posteriori_backward_error_testX =  norm(b-A*testX,1) / (norm(A,1) * norm(testX,1) * norm(b,1));
    exact_forward_error_testX = norm(testX - x_sol,1) / norm(x_sol,1);
    fprintf('----Sfalmata gia to sistima x = A \\ b, gia to mitroo %s-----', mx_id);
    fprintf('\nFragma Dikti Katastasis: %d',cond_A);
    fprintf('\nSinepagwmeno Empros Sfalma: %d\n',sinepagomeno_forward_error_testX);
 	fprintf('A Posteriori Sxetiko Pisw Sfalma: %d\n',posteriori_backward_error_testX);
 	fprintf('Akrives sxetiko Pisw Sfalma: %d\n', exact_forward_error_testX);
    
    %----------O PIO PANW KWDIKAS ITAN MEXRI KAI TO ERWTIMA 3---------
    %-----------------------------------------------------------------
    %-----------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %----------TIS XRONOMETRISEIS KAI TWN ELEGXO APOTELESMATWN TA EXW VALEI
    %----------SE COMMENTS, AFOU EKTELESTIKAN MONO GIA TESTING-------------
    
    
    
    %--------------------------------ERWTIMA 4-----------------------------
    blk = 16;
    X_block = SMW_solve_0052_block(A,b,M,P,Q,blk)
    
    %ERWTIMA 4.b
    empros_sfalma_block = norm(X_block - x_sol,1);
    fprintf('\nEmpros Sfalma (periptwsi block): %d\n',empros_sfalma_block); 
    

    %Xronometrisi algorithmou SMW_solve_0052_block gia ta
    %diafora blk (erwtima 4.c)
    %blk=1;
    %blk=2;
    %blk=4;
    %blk=8;
%     blk=16;
%     for i=1:100
%         tic
%         X = SMW_solve_0052_block(A,b,M,P,Q,blk);
%         t1(i)=toc;
%     end
%     avg_xronos_ektelesis = sum(t1)/100;
%     fprintf('Mesos Xronos ektelesis gia blk=%d: %f\n',blk, avg_xronos_ektelesis);
     
    %Xronometrisi arxikou algorithmou SMW_solve_0052
    %xriastike sto erwtima 4.c
%     for i=1:100
%         tic
%         X = SMW_solve_0052(A,b,M,P,Q,sdir);
%         t1(i)=toc;
%     end
%     avg_xronos_ektelesis_arxikou = sum(t1)/100;
%     fprintf('Mesos Xronos ektelesis arxikou algorithmou: %f\n', avg_xronos_ektelesis_arxikou);




    %--------------------------------ERWTIMA 5-----------------------------
    b=b';
    X_mex = SMW_solve_0052_block_mex(A,b,M,P,Q,blk)


    %Elegxos swstwn apotelesmatwn apo tin sinartisi
    %SMW_solve_0052_block_mex (erwtima 5.a)
%     for i=1:100
%         [X]=SMW_solve_0052_block(A,b,M,P,Q,blk);
%         b=b';
%         [X1]=SMW_solve_0052_block_mex(A,b,M,P,Q,blk);
%         b=b';
%         nX(i)=norm(X-X1);
%     end
%     norm(nX)
    
    %Xronometrisi algorithmou SMW_solve_0052_block_mex gia ta
    %diafora blk (erwtima 5.b)
%     b=b'
%     for i=1:100
%         tic
%         [X] = SMW_solve_0052_block_mex(A,b,M,P,Q,blk);
%         t1(i)=toc;
%     end
%     avg_xronos_ektelesis_mex = sum(t1)/100;
%     fprintf('Mesos Xronos ektelesis MEX gia blk=%d: %f\n',blk, avg_xronos_ektelesis_mex);

end