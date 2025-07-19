# Virtual-Flagellum
Code to simulate a free-swimming Internally-driven Kirchhoff Rod (IDKR) in bulk fluid and near a plane wall. Matlab code suitable for MEX-function generation. 
The main functions are main_RPY and main_Ani. To run, please open Example_Script.m.

Both main_RPY and main_Ani are suitable for MEX-function generation using Matlab Coder. 
The following two lines generate a MEX-function for main_RPY called main_RPY_mex, which can be called in place of main_RPY. Please ensure you have Matlab Coder available.
Before running the lines below, please comment out lines 123-137 in main_RPY (or main_Ani) as plotting is not supported during MEX-function generation.

```
ins = {0.01,0.02,50,[0 0 9],[0 0 10],[0 0 4*pi],1,1,1,0,1e-8,0.01,10,1e-6} % Example inputs to main_RPY
codegen main_RPY -args ins
```

Following this, in Example_Script.m, the function call "main_RPY" can be replaced with "main_RPY_mex" to achieve a much shorter running time.
