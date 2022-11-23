In order to run these files (in particular realSheepSMC.R found in the folder Rcode)
then you need to download the necessary packages
(e.g. install.packages("mcmcse") )
and need to run in the terminal:
sudo apt-get install libgsl-dev
but also, you need to compile the dissPackage3 c++ functions using:
'cd /folder-location/MoutonModel/'
where the subfolders should be 'dissPackage3', 'Plots', 'Rcode', etc. Then do
'R CMD build dissPackage3'
then
'R CMD INSTALL dissPackage3'
This will insert the package into the user executables for Linux OS.

Additionally, run:
sudo apt-get install libgsl-dev

















For GPU programming in R such as through gmatrix, using HWP GPU updates to SMC code.
If you are lucky, this might work: 'https://www.rdocumentation.org/packages/gmatrix/versions/0.3'
   ......That probably didn't work....... so here goes:

1) PRAY TO GOD: IF YOU DON'T BELIEVE, NOW IS THE TIME TO HAVE HER/HIM ON YOUR SIDE

2) Back-up everything on your computer, just in case... SERIOUSLY, DO IT!

3) Ensure you have an Nvidia graphics card: 'lspci | grep -i nvidia'
   IF NOT: YOU SHALL NOT PASS! Change computer

4) Update your Nvidia graphics driver:
   	  'sudo apt update && sudo apt upgrade'
   	  'sudo apt remove nvidia*'
   	  'sudo ubuntu-drivers autoinstall'

5) Find out the driver currently installed: 'nvidia-smi'
   EXAMPLE: I have driver 435.21 and CUDA 10.1
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 435.21       Driver Version: 435.21       CUDA Version: 10.1     |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|===============================+======================+======================|
|   0  GeForce 940M        Off  | 00000000:04:00.0 Off |                  N/A |
| N/A   44C    P0    N/A /  N/A |      0MiB /   983MiB |      1%      Default |
+-------------------------------+----------------------+----------------------+
+-----------------------------------------------------------------------------+
| Processes:                                                       GPU Memory |
|  GPU       PID   Type   Process name                             Usage      |
|=============================================================================|
|  No running processes found                                                 |
+-----------------------------------------------------------------------------+
	### IN CASE IT ERRORS: 'sudo prime-select nvidia'
	then retry... you have now set the NVIDIA GPU in use instead of Intel on-board GPU

6) If CUDA is not installed but the drivers are, use: 'sudo apt install nvidia-cuda-toolkit'
   It may be wise to check that the driver and CUDA version are compatible using table 1 of
   https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html

7) For information on your specific GPU, including which sm_50 go to:
   https://developer.nvidia.com/cuda-gpus
   Where the number in sm_50 is because I have the NVIDIA GeForce 940M
   which has a compute capacity of 5.0 -> sm_50

8) Finally, download gmatrix from the github using the following:
   i)   'cd ~/Downloads'
   ii)  'git clone https://github.com/njm18/gmatrix.git'
   iii) 'rm -rf ./gmatrix/.git'
   iv)  'MAKE="make -j7"'
   v) 	'CUDA_HOME="/usr/local/cuda"'
   vi)	'R CMD INSTALL gmatrix --configure-args="--with-arch=sm_50 --with-cuda-home=$CUDA_HOME"'

9) You should now be able to access this in Rstudio using 'library("gmatrix")'

10) BE CAREFUL USING 'sudo prime-select nvidia': THIS CAN MESS YOUR GRAPHICS UP
    (especially for a sexy multiple screen setup)...
    Use 'prime-select query' to check that NVIDIA is only used when necessary.
    
Good luck, my young padawan. May the GeForce be with you (I'm sorry, I'm so sorry).
   
