function [Filament, t] = main_Ani(a,ds,Ns,SN,m_0,k_b,phi,mu,KB,KT,wall,dtfac,num_settling_times,outfilename,save_step,concheck_tol) %#codegen
%   Supplementary code to 'Elastohydrodynamic mechanisms govern beat pattern transitions in eukaryotic flagella', by S Veeraragavan, F Yazdan Parast, R Nosrati, and R Prabhakar. 
%   Access the paper at https://doi.org/10.1101/2024.02.04.578806.
% 
%   This code is built over the numerical methods developed and published by Schoeller et al. (J. Comp. Phys. 424, 2021), available at https://doi.org/10.1016/j.jcp.2020.109846. 
%
%   It simulates the three-dimensional motion of a single free-swimming flagellum (modelled as an internally-driven Kirchhoff rod) immersed in viscous fluid.

% Toggles for plotting and printing results to the terminal   
print_results = true
plot_filament = true

% Filament data
weight_per_unit_length = 0;         % weight per unit length of the filament 
L = ds*Ns;                          % total length            
omega_b = [SN(1)^4/L^4/mu*KB, SN(2)^4/L^4/mu*KB, SN(3)^4/L^4/mu*KB];       % Beat frequency (3 components in the internal frame: {d1, d2, d3})

% Initialise Filament
Filament = struct('N_w',Ns,'KB',KB,'KT',KT*KB,'DL',ds,'Length',L,'StrainTwist',zeros(3,Ns-1),'weight_per_unit_length',weight_per_unit_length,'R',a*ones(1,Ns),'X',zeros(3,Ns),'Xm1',zeros(3,Ns),'Xm2',zeros(3,Ns),'X1',zeros(6,1),'Q',zeros(Ns,12),'U',zeros(9,Ns),'V',zeros(6,Ns),'F',zeros(6,Ns),'Lambda',zeros(3,Ns-1),'Lmat',zeros(6*Ns,6*Ns),'Umat',zeros(6*Ns,6*Ns));
Filament = InitialSetup(Filament,[0;0;wall]);

% Choose timescale
unit_time = mu*L^4/KB; % 1 unit time: viscoelastic timescale
dtve = mu*(ds)^4/KB;   % Smallest viscoelastic timescale
dtva = mu*(ds)^3*L/KB/max(m_0); % Smallest viscous-active timescale
dt = min([dtve dtva]);
steps_per_unit_time = unit_time/dt;
TOTAL_STEPS = floor(num_settling_times*steps_per_unit_time);

t = 0;
save_now = save_step - 1;

% Time and iteration counts
N = 6*Filament.N_w;
N = [0,cumsum(N)];
Nbroy = N(end);
max_broyden_steps = 5*Ns;
frame_time = zeros(TOTAL_STEPS,1);
iters = zeros(TOTAL_STEPS,1);       % Number of Broyden's iterations
running_total_count = 0;            % For average number of Broyden's iters

% Time Integration - Main Simulation Loop

i = 1;
for nt = 1:TOTAL_STEPS  
    iter = 0;

    p_broy = max_broyden_steps + 1;
    Cmat = zeros(Nbroy,p_broy); 
    Dmat = zeros(Nbroy,p_broy); 
    
    frame_start = tic;    
    
    % Update positions based on initial guess
    Filament = InitialGuess(Filament);
    Filament = RobotArm(Filament);
   
    % Calculate the active moment vector
    Filament = ActiveST(Filament,omega_b, k_b, phi, m_0, t);
    
    % Find f(X_k) and place into ERROR_VECk, the error vector
    % If ||ERROR_VECk|| < concheck_tol, then concheck = 0. Else, 1.    
    [concheck,ERROR_VECk,Filament,FS] = F_Ani(Filament,mu,dt,nt,concheck_tol,wall); 
    
    % Update the guess using the Bad Broyden method, until ||ERROR_VECk|| < concheck_tol (concheck = 0)
    % See paper by Schoeller et al. (referenced in header) for details.
    
    Filament = ConstructAndDecomposeJacobian(Filament,dt,mu);
    J0invERROR_VECk = blockwise_backslash_jacobian(Filament, ERROR_VECk);    
     
    num_broydens_steps_required = 0;
    while (concheck == 1) % Alg 2, Line 4
        % Alg 2, Line 5. DeltaX is Delta X in paper.
        DeltaX = -apply_inverse_jacobian(J0invERROR_VECk, Cmat, Dmat,...
                                         ERROR_VECk, iter);

        % Update the positions and lambdas
        Filament = ApplyUpdate(Filament,DeltaX(N(i)+1 : N(i+1)));
        Filament = RobotArm(Filament);       
        
        % Check to see if the new state is an acceptable solution:
        % ERROR_VECk1 = f(X_(k+1))        
        [concheck,ERROR_VECk1,Filament,FS] = F_Ani(Filament,mu,dt,nt,concheck_tol,wall);     
               
        iter = iter + 1;
        
        y_vec = ERROR_VECk1 - ERROR_VECk;
        J0invERROR_VECk1 = blockwise_backslash_jacobian(Filament, ...
                                                              ERROR_VECk1);
        
        y_vec_sq = y_vec'*y_vec;
        Cmat(:,iter) = -apply_inverse_jacobian(J0invERROR_VECk1, ...
                                            Cmat, Dmat, ERROR_VECk1, iter);
        Dmat(:,iter) = y_vec/y_vec_sq;
        ERROR_VECk = ERROR_VECk1;
        J0invERROR_VECk = J0invERROR_VECk1;
        
        % If the number of iterations maxes out, proceed to next timestep
        % anyway and see what happens (but flag it with a *)
        if (iter > max_broyden_steps)
            fprintf(' *');
            concheck = 0;
        end

        num_broydens_steps_required = num_broydens_steps_required + 1;
        running_total_count = running_total_count + 1;        
    end
    
    % Step forward in time
    t = t + dt;    
    Filament = EndOfStepUpdate(Filament);
    
    % Print results to terminal
    save_now = save_now + 1;    
    if(save_now == save_step && print_results)
        if nt == 1 
          fprintf('dt, a, Ns, L, m0, SN, k_a, phi, mu, KB, KT, dtfac, tol \n');
          fprintf('%.4e, %.6f, %.0f, %.4f, %.0f, %.0f, %.0f, %.0f, %.2f, %.2f, %.2f, %.4f, %.4f, %.4f, %.0f, %.0f, %.0f, %.0f, %.0f, %.0f, %.0f, %.2e \n',dt, a, Ns, L, m_0(1), m_0(2), m_0(3), m_0(4), SN(1), SN(2), SN(3), k_b(1), k_b(2), k_b(3), phi(1), phi(2), phi(3), mu, KB, KT, dtfac, concheck_tol);
        end
        PrintToFile(Filament,t);
    end    
    
    if save_now == save_step
        save_now = 0;
    end

    % Plot
    if plot_filament 
       com = mean(Filament.X,2);
       L = Filament.Length;
       x = Filament.X;
       plot3(x(1,:)/L,x(2,:)/L,x(3,:)/L,'-','LineWidth',10);
       pbaspect([1 1 1])
       xlim([com(1)/L-0.5,com(1)/L+0.5]);
       ylim([com(2)/L-0.5,com(2)/L+0.5]); 
       zlim([com(3)/L-0.5,com(3)/L+0.5]); 
       xlabel('(x-x_{COM})/L');
       ylabel('(y-y_{COM})/L');
       zlabel('(z-z_{COM})/L');
       pause(0.001);
    end
    
    frame_time(nt) = toc(frame_start);
    iters(nt) = iter;
end

fprintf('\n');
disp('Run finished');
disp(['Total time:' format_time(sum(frame_time))]);
end
