function out = InvertLocalBlock(Fil,v)
    % Filament.INVERTLOCALBLOCK(v)  returns J_0\v, where J_0 is the
    %                               pre-computed, pre-LU-decomposed
    %                               approximate Jacobian for that 
    %                               filament.
        out = Fil.Lmat\v;
        out = Fil.Umat\out;
    end