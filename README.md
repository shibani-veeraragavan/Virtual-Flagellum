# Virtual-Flagellum
Code to simulate a free-swimming Internally-driven Kirchhoff Rod (IDKR) in bulk fluid and near a plane wall. Matlab code suitable for MEX-function generation. 
The main functions are main_RPY and main_Ani. To run, please open Example_Script.m.

Both main_RPY and main_Ani are suitable for MEX-function generation using Matlab Coder. 
The following two lines generate a MEX-function for main_RPY called main_RPY_mex, which can be called in place of main_RPY. Please ensure you have Matlab Coder available.
Before running the lines below, please comment out lines in main_RPY (or main_Ani) as plotting is not supported during MEX-function generation.

  ins = {a,ds,Ns,S,A,k,mu,KB,KT,wall,dt,0.01,save_step,concheck_tol}
  
  codegen main_RPY -args ins

Following this, in Example_Script.m, line, "main_RPY" can be replaced with "main_RPY_mex" to achieve a much shorter running time.
