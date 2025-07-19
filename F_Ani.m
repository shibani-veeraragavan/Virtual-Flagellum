function [concheck_local,ERROR_VECk1_local,Filament] = F_Ani(Filament,mu,dt,nt,tol,wall)
% F  places forces and torques on the segments, calculates the resultant
%    velocities and angular velocities, and forms the error vector f(X*).
%    Then checks convergence. For details, see docstrings of functions
%    within.

    Filament = InternalForcesAndTorques(Filament,dt);
    
    [Filament, FS] = collision_barrier(Filament);
    
    Filament = Ani_Wall(Filament,mu,wall);
    
    % Check convergence between x_(n+1) and x_n, and also check the
    % constraint. concheck = 0 if all fine, 1 otherwise. The error vectors
    % are all compiled into ERROR_VECk1_local.    
    [concheck_local,ERROR_VECk1_local] = constraint_check_robot_arm(...
                                                      Filament,dt,nt,tol);
end
