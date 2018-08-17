#!/bin/python python
#Get velocity of a droplet as a function of time, from 'film_xmol' file

def old_main(): # main_script
    """
    This routine reads a film_xmol, and writes the position of the droplet as a function of time
    """
    import numpy as np
    import sys

    if len(sys.argv) == 2: # Check arguments passed when the script is called
        Dx_win = float(sys.argv[1])  # Define  length in x to obtain X_cm for the droplet. The narrower the window
                # the better resolution is obtained, given that it's wider than droplet's diameter
    else:
        print "ERROR: Script must be called with one argument."
        print "The (aproximate) largest diameter of the dropplet must be passed"
        print "Example: $ python get_cmPos_vs_t.py 40. "     
        print ""
        quit()

    Lx, Ly, Lz, frame_tot, n_part, n_b, dt =  read_params() # read system paramets 
    
    f = open('refined_cm_vs_t.mide', 'w+')

    for frame_nr in range(1,int(frame_tot+1)):   # Loop over frames of film_xmol
        r0  =  read_conf2(frame_nr, n_part, n_b) # Read liquid configuration, for this frame

        cm_mode_x = histo_mode(r0[0,:],Lx, n_part) # Get first cm guess extracting the mode 
        mode=np.copy(cm_mode_x)                                          # of the histogram of melt positions in x
       
        if Lx < Dx_win: # IF Simulation Box is smaler than Ds_win, get droplet's mass center directly  
            pos_cm = center_drop(r0,Lx, n_part, cm_mode_x)
             
        # If Lx is large, do a refined calculation of the droplet's cm position 
        elif Lx > Dx_win :  # Inaccuracy increaces due to large gas/liquid ratio of particles
            pos_cm = np.zeros(3)
            pos_cm[0] = refined_pos_cm_x(r0[0,:], Lx, cm_mode_x, Dx_win)
            pos_cm[1], pos_cm[2] = refined_pos_cm_yz(r0, Lx, pos_cm[0], Dx_win)

        # Save step of trajectory to file f : time, x_cm, y_cm, z_cm
        s = str(str(frame_nr*dt)+' '+str(pos_cm[0])+' '+str(pos_cm[1])+' '+str(pos_cm[2])+ '\n')     
        f.write(s) 

    f.close()

def read_conf2(frame_nr, n_liq, n_brush):
    """This routine reads the positions of the melt and brush for a specific frame from file 'film_xmol'.
    The outputs are two numpy arrays x0 (melt) and x0_b(brush) with the dimension in the first index, and particle number
    in the second x0[0][0] is the x position of the first particle.
    This function should be calles as: r0 = read_conf(1); 
    to obtain the first frame (frame 0 does not exist)"""
    import numpy as np
    import linecache
    x0 = np.empty([3,n_liq]) # defie a numpy array of size 3 x melt_pos

    _split = str.split 
    frame_lines = int(n_liq + n_brush + 2) # Select lines to read from film_xmol
    line_init_m = int((frame_nr - 1) * frame_lines + 3 + n_brush) 
    line_end_m  = int(line_init_m + n_liq)

    j=0
    for i in range(line_init_m,line_end_m):# Save melt positions in another array
        x = _split(linecache.getline('film_xmol', i))[1:4]
        x = np.array(x, dtype='|S16')
        x = x.astype(np.float)
        x0[:,j] = x[:]
        j+=1
        
    return x0 # return are two  numpy arrays with melt particles in 'em   

def center_drop(r0, Lx, n_part, x_cm_guess):
    """
    This routine takes a configuration in r0, and a box length Lx.
    It calculates the center of mass of a droplet. It does not work for more than one droplet.
    The idea is to calculate the cm of the system; then transalte the system so that the cm is in the middle of
    the symulation box; apply PBC; repeat. When the algorithm congerges, it is possible to know the center of mass
    of the droplet, by summing the translatios.
    The first argument must be a numpy array with the positions
    The second argument must be the box length
    3rd argument: # of particles
    4th argument: The first guess of x_cm of droplet. Works best with mode of particle position histogram. Otherwise use 
    system center of mass
    """
    import numpy as np
    r0_trans = np.copy(r0) #FUCKING IMPORTANT TO DO A COPY, OTHERWISE IT OVERWRITES r0
    pos_cm = r0_trans.sum(1) / float(n_part)
    x_cm = x_cm_guess    #OLD x_cm = (r0_trans.sum(1) / r0_trans.shape[1])[0]

    x_trans = 0 # Translation amount in x

    while True: # Iteratice proces to center droplet in the box
        delta_x = x_cm - Lx/2.
        r0_trans[0,:] = r0_trans[0,:] - delta_x # Carry out translation         
        r0_trans[0,:] = r0_trans[0,:] - Lx*( 2. * r0_trans[0,:] / Lx - 1).astype(int)    #Implement PBC in x
        x_cm_new = (r0_trans.sum(1) / r0_trans.shape[1])[0] #Calculate new centre of mass
        delta_x_new = x_cm_new - Lx/2.

        if np.abs(delta_x) < 0.1: # Convergence condition 1
            if np.abs(delta_x_new) <= np.abs(delta_x): # Convergence condition 2
                break

        x_trans = x_trans +  delta_x # Accumulate translations done
        x_cm = x_cm_new # Remember old center of mass position for next iteration

    pos_cm[0] = x_trans + Lx/2. # Calculate initail position of droplet
    pos_cm[0] = pos_cm[0] - Lx * ( 2. * pos_cm[0] / Lx - 1. ).astype(int) # Apply PBC
    return  pos_cm


def center_drop_x(x0,Lx,n_part):
    """
    This routine rakes a the x components of a configuration in r0, and a box length Lx.
    It translates the centre of mass to x = Lx/2
    The first argument must be a numpy array with the positions only in x
    The second argument must be the box length
    The third argument, the number of particles
    """
    import numpy as np
    r0_trans = np.copy(x0) # PYTHON SUBTLETY: r0_trans AND x0 ARE NOT SUPPOSED TO POINT TO THE SAME MEMORY SPACE
                           # THEY ARE SEPARATE VARIABLES
    x_cm = r0_trans.sum(0) / float(n_part)
   
    x_trans = 0 # SEGUIR DEFINIENDO TANS FINAL

    while True:
        delta_x = x_cm - Lx/2.
        r0_trans[:] = r0_trans[:] - delta_x # Carry out translation
        #r0_trans[0,:] = apply_PBC(r0_trans[0,:],Lx,n_part) # Implement PBC slow 
        r0_trans[:] = r0_trans[:] - Lx*( 2. * r0_trans[:] / Lx - 1).astype(int)    #Implement PBC in x
        x_cm_new = r0_trans.sum(0) / float(n_part) #Calculate new centre of mass
        delta_x_new = x_cm_new - Lx/2.

        if np.abs(delta_x) < 0.1:
            if np.abs(delta_x_new) <= np.abs(delta_x):
                break

        x_trans = x_trans +  delta_x
        x_cm = x_cm_new

    pos_cm = x_trans + Lx/2.
    pos_cm = pos_cm - Lx * int( 2. * pos_cm / Lx - 1. )
    return  pos_cm
    
def read_params():
    """
    This script reads Lx, Ly, Lz, n_frames, nm, n_brush, delta_t from system_input and mfa_input
    """
    import numpy as np
    import os.path

    for file_name in 'mfa_input', 'system_input':
        if os.path.isfile(file_name) == False:
            print "Fatal Error: File "+ file_name +" is missing"
            print "Finishing script here"
            exit()

    j=0
    with open('mfa_input') as f:
        for line in f:
            j += 1
            if(j==3):
                n_steps = float(line.split()[0])
            if(j==4):
                n_obs = float(line.split()[0])    
            if(j==5):
                dt = float(line.split()[0])                                     
            if(j==12):
                Lx = float(line.split()[0])
            if(j==13):
                Ly = float(line.split()[0])
            if(j==14):
                Lz = float(line.split()[0])
    delta_t = dt * float(n_obs)    
    n_frames = int(n_steps / n_obs)                      



    #Now read system_input
    j=0
    with open('system_input') as f:
        for line in f:
            j += 1
            if(j==5):
                Lx = Lx * float(line.split()[2])
                Ly = Ly * float(line.split()[3])
                n_mon = float(line.split()[0])
                n_chain = float(line.split()[1])
            if(j==7):
                n_mon_d = float(line.split()[0])
                n_chain_d = float(line.split()[1])



    #j=0
    #with open('conf_old') as f:
    #    for line in f:
    #        j += 1
    #        if(j==3):
    #            Lx = Lx * float(line.split()[0])
    #            Ly = Ly * float(line.split()[1])

    #        if(j==2):    
    #            n_mon = float(line.split()[0])
    #            n_chain = float(line.split()[1])
    #            n_mon_d = float(line.split()[2])
    #            n_chain_d = float(line.split()[3])

    return Lx,Ly,Lz,int(n_frames),int(n_chain_d*n_mon_d),int(n_mon*n_chain), delta_t


def read_conf_brush(frame_nr, n_liq, n_brush, run_nr):
    """This routine reads the positions of the melt and brush for a specific frame from file 'film_xmol'.
    The outputs are two numpy arrays x0 (melt) and x0_b(brush) with the dimension in the first index, and particle number
    in the second x0[0][0] is the x position of the first particle.
    This function should be calles as: r0 = read_conf(1); 
    to obtain the first frame (frame 0 does not exist)"""
    import numpy as np
    import linecache
    import os
    import gzip

    file_tgt = str(str(run_nr)+"_run/"+"film_xmol")

    if os.path.isfile(file_tgt) == False:

        if os.path.isfile(file_tgt+".gz") == False:
            print "Fatal Error: File "+ file_name +" is missing"
            print "Finishing script here"
            exit()

        else: # unzip film_xmol            
            inF = gzip.GzipFile(file_tgt+".gz", 'rb')
            s = inF.read()
            inF.close()
            
            outF = file(file_tgt, 'wb')
            outF.write(s)
            outF.close()

            os.remove(file_tgt+".gz")

    x0 = np.empty([3,n_brush]) # defie a numpy array of size 3 x melt_pos

    _split = str.split 
    frame_lines = int(n_liq + n_brush + 2) # Select lines to read from film_xmol
    #line_init_m = int((frame_nr - 1) * frame_lines + 3 + n_brush) 
    line_init_b = int((frame_nr - 1) * frame_lines + 3 ) 
    
    #line_end_m  = int(line_init_m + n_liq)
    line_end_b  = int(line_init_b + n_brush)

    j=0
    for i in range(line_init_b,line_end_b):# Save melt positions in another array
        x = _split(linecache.getline(file_tgt, i))[1:4]
        x = np.array(x, dtype='|S16')
        x = x.astype(np.float)
        x0[:,j] = x[:]
        j+=1

    #print "first line read from "+file_tgt+" is ",line_init_b
    #print "last line read from "+file_tgt+" is ",line_end_b-1    

    return x0 # return are two  numpy arrays with melt particles in 'em   


 
def main():
    """ get 2D g(r) for the brush particles """
    import numpy as np
    from gr_lib import gr_engine as calc_gr
    
    Lx, Ly, Lz, n_frames, n_liq, n_brush, delta_t = read_params()

    # Set default parameters
    delta_z = 0.5   # height of slab to take 2D g(r)
    r_bin = 0.1     # resolution in the radial direction

    # Get first and last run directories, to analyze film_xmol
    first_run, last_run = get_run_dirs()
    
    # Set parameters according to what is selected above
    z1 = Lz/2. - delta_z / 2.0
    z2 = Lz/2. + delta_z / 2.0 
    n_bins = int( .5 * np.sqrt( np.square(Lx) + np.square(Ly) ) / r_bin ) + 1
    g_r = np.zeros(n_bins)    
   
    # Loop over all the simulations and get statistics for g(r):
    n_fr_tot = 0
    for run_nr in range( first_run, last_run + 1 ):
        for i_frame in range(1, n_frames + 1 ):
            n_fr_tot += 1
            r0 = read_conf_brush(i_frame, n_liq, n_brush, run_nr)
            g_r_tmp = calc_gr(n_bins, r_bin, z1, z2, Lx, Ly,r0, n_brush) # Fortran enginge is summoned
            g_r += g_r_tmp

    # Get mean value of g(r)
    g_r = g_r / float(n_fr_tot)
    
    save_results("2D_gr.mide",g_r, n_bins, r_bin)

def save_results(name, vector, dim, delta_x):

    f = open(name, 'w+')

    for i in range (dim):
        s = str( str( delta_x * float(i) ) + " " + str(vector[i]) + '\n')
        f.write(s)

    f.close()

def get_run_dirs():
    import os
    last_run = 1
    first_run = 999999999
    for dirs in os.listdir("./"):
        if os.path.isdir(dirs) == True:
            if dirs.endswith("_run") == True:
                run_nr = int( dirs.split("_")[0] )
                if run_nr > last_run:
                    last_run = run_nr
                if run_nr < first_run:
                    first_run = run_nr
                
    if first_run < 20 and last_run > 20:
        first_run = 20
 
    return first_run, last_run

def refined_pos_cm_x(x0, Lx, x_cm_guess, Dx_win):
    """ Do refined droplet cm position calculation. It's the same as center_droplet, but cutting
    the length of the box in iterative steps, to reduce the gas_particles / liquid_particles ratio.
    This should increase the accuracy of the meassurments"""
    import numpy as np

    Delta_x_2 = np.copy(Dx_win)  # Define half window in X for refined meassurment of X_cm
    if 2*Delta_x_2 > Lx:
        Delta_x_2 = Lx/2. 
    #Loop here to compact code      
        
    #Get particles that  x_cm_guess - Delta_x_2 < x < x_cm_guess + Delta_x_2 i
    r0_aux = get_win_part(x0, x_cm_guess, Lx, Delta_x_2)               
    part_in_win = r0_aux.shape[0]   # how many particles are in refined window              
    
    # Get x_cm refined for this window in x   
    x_correction = center_drop_x(r0_aux, 2.*Delta_x_2, part_in_win) #Get correction in x_cm
    x_cm_guess = x_cm_guess + x_correction - Delta_x_2
    x_cm_guess = x_cm_guess - Lx * int( 2. * x_cm_guess / Lx - 1.) # Apply PBC 
         
    Delta_x_2 = Dx_win/2. # Narrow Refined Window
    
    # Make final refined meassurment of x_xm_guess, for a window of Dx_win
    r0_aux = get_win_part(r0_aux, x_correction, 2.*Delta_x_2, Delta_x_2)
    part_in_win = r0_aux.shape[0]
    # Get x_cm refined for this window in x   
    x_correction = center_drop_x(r0_aux,  Dx_win , part_in_win) #Get correction in x_cm    
    x_cm_guess = x_cm_guess + x_correction - Delta_x_2
    x_cm_guess = x_cm_guess - Lx * int( 2. * x_cm_guess / Lx - 1. ) # Apply PBC          
        
    return  x_cm_guess

def refined_pos_cm_yz(r0 ,Lx, x_cm, Dx_win):
    """ Get cm positions in y and z, for refined window
    """
    import numpy as np

    Delta_x_2 = Dx_win / 2.  
    y_cm = 0.
    z_cm = 0.
    n_part_in_win = 0 
    i_part = 0
    
    for x_part in r0[0,:]: #Loop over all x_coorinates of particles. 
        dis_x = x_part - x_cm - Lx * int( 2. * ( x_part - x_cm ) / Lx) # calc dis in x
        
        if (dis_x < Delta_x_2 and dis_x > - Delta_x_2): # IF particle is inside refined window
            y_cm += r0[1,i_part]
            z_cm += r0[2,i_part]
            n_part_in_win += 1                    
        
        i_part += 1
    
    y_cm = y_cm / float( n_part_in_win ) 
    z_cm = z_cm / float( n_part_in_win )
    
    return y_cm, z_cm

def get_win_part(x0, x_cm_guess, Lx, Delta_x_2):
    """ Returns a numpy array with all ther particles of x0 inside a box of length 2 * Delta_x_2, and
    centered in x_cm_guess. 
    The new vector starts at x = 0 and ends in x= 2 * Delta_x_2
    """
    import numpy as np

    r0_aux_ls = [] 
 
    for x_part in x0: #Loop in all particles. 
        dis_x = x_part - x_cm_guess - Lx * int( 2. * ( x_part - x_cm_guess ) / Lx) # calc dis in x
        
        if (dis_x < Delta_x_2 and dis_x > - Delta_x_2): # IF particle is inside refined window
            r0_aux_ls.append(x_part)                
   
    r0_aux = np.asarray(r0_aux_ls) - x_cm_guess + Delta_x_2 #r0 as numpy array and translation
                    
    i_part = 0 
    for x_part in r0_aux:  #Correct for PBC in Lx, so that all particles lie in little window between 0 and Delta_x 
        if x_part  < 0.: 
            r0_aux[i_part] = x_part + Lx 
        elif r0_aux[i_part] >  Lx:
            r0_aux[i_part] = x_part - Lx 
        i_part += 1
    
    return r0_aux        

    
def histo_mode(x0,Lx, n_part):
    """ Do a histogram of particles positions in x, and return the mode
    """
    import numpy as np

    bin_w = 1.0 # Define width of x bins
    n_bins = int( Lx / bin_w ) + 1 # bin number
    histo_x = np.zeros(n_bins, dtype=int) # Define histogram
    index = ( x0 / bin_w ).astype(int) # calculate indexes for histogram

    for i in index: # Fill histogram
        histo_x[i]  += 1
    
    x_mode = np.argmax(histo_x) # Get mode from histogram
    return  x_mode # return mode from distribution


########### END SUPPORT ROUTINES ##########

if __name__ == "__main__":
    main()

