function x = SMW_solve_0052(A,b,M,P,Q,sdir)
    m = size(A);
    n = m(1); %dimension of A
    x = M \ b; %find x, by solving Mx=b
    y = M \ P; %find y, by solving My=P
  
    for l=1:n-1
    	x = x-((Q(:,l)' * x) / (1+Q(:,l)' * y(:,l)))*y(:,l); %compute x
        for k=(l+1):n
            y(:,k) = y(:,k)-(Q(:,l)'*y(:,k)/(1+Q(:,l)'*y(:,l)))*y(:,l);
        end
    end
        x = x-(Q(:,n)'*x/(1 + Q(:,n)' * y(:,n)))*y(:,n); %compute final result of vector x                     
end