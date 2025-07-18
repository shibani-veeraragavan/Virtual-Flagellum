function PrintToFile(Fil,t)
    % Filament.PRINTTOFILE(fid,t)  prints a string of filament segments 
    %                              positions, quaternions, velocities and 
    %                              forces to the terminal.

        fprintf('\n');
        for j_w = 1:Fil.N_w
            if j_w == Fil.N_w
                L1 = 0;
                L2 = 0;
                L3 = 0;
            else
                L1 = Fil.Lambda(1,j_w);
                L2 = Fil.Lambda(2,j_w);
                L3 = Fil.Lambda(3,j_w);                    
            end
            
            
            fprintf(['%.8f' repmat(', %.6f',1,22)], ...
            t, Fil.X(1,j_w), Fil.X(2,j_w), Fil.X(3,j_w), ...
            Fil.Q(j_w,1), Fil.Q(j_w,2), Fil.Q(j_w,3), Fil.Q(j_w,4), ...
            Fil.V(1,j_w), Fil.V(2,j_w), Fil.V(3,j_w), ...
            Fil.V(4,j_w), Fil.V(5,j_w), Fil.V(6,j_w), ...
            Fil.F(1,j_w), Fil.F(2,j_w), Fil.F(3,j_w), ...
            Fil.F(4,j_w), Fil.F(5,j_w), Fil.F(6,j_w), ...
            L1, L2, L3);
        fprintf('\n');
        end
end
