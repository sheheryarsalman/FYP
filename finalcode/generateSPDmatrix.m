function A = generateSPDmatrix(n,covSF)
% Generate a dense n x n symmetric, positive definite matrix

    A = covSF*rand(n,n); % generate a random n x n matrix taking val.
                                  % from 0 to covScaleFactor
                                  
    A = 0.5*(A+A');               % construct a symmetric matrix
    
                                  % add random values to ensure PD
    A = A + diag(diag(randi([10 10*covSF], n)/(10/covSF)));                                   
    [~,p] = chol(A);              %test for positive-definiteness
    if(p!= 0)
        A = A + diag(diag(randi([10 25], n)/(10/covSF))); 
        warning('Matrix may not be PD. Consider running again.')
    end
end      