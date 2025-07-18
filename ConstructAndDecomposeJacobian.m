 function Fil = ConstructAndDecomposeJacobian(Fil,dt,mu)
    % Filament.CONSTRUCTANDDECOMPOSEJACOBIAN(dt,mu)  constructs the
    %     approximate Jacobian for that filament. It places the
    %     decomposed parts into Filament.Lmat and Filament.Umat.
    approxj = approximate_jacobian(Fil,dt,mu);
    [L,U] = lu(approxj); 
    Fil.Lmat(:,:) = L;
    Fil.Umat(:,:) = U;
    end