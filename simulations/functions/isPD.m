function tf = isPD(A)
    % Symmetrize first
    A = (A + A')/2;      
    [~,p] = chol(A);
    tf = (p == 0);        % true if PD
end