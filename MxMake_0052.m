function A = MxMake_0052(mx_id,n)
    switch mx_id
        case 'had'   
            A=hadamard(n);
        case 'trihad'
            A=triu(hadamard(n));
        case 'toep'
            A=toeplitz([4,-1,zeros(1,n-2)])
        case 'mc'
            s=2; th=1;
            for i=1:n
                A(i,i) = 1 + i.^th;
                for j=1:n
                    if (i~=j)
                        A(i,j)=1/(abs(i-j)^s);
                    end
                end 
            end
        case 'wathen'
            A=gallery('wathen',n,11);
        case 'CollegeMsg'
            struct_CollegeMsg = load('CollegeMsg.mat');
            A=struct_CollegeMsg.Problem.A;
            A=eye(1899)-0.85*A;
    end
end