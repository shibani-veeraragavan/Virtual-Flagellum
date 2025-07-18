function com = CentreOfMass(Fil)
    % Filament.CENTREOFMASS()  returns the centre of mass of the
    %                          filament.
        com = mean(Fil.X,2);
    end