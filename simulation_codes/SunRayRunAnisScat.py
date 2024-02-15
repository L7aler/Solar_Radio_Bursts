# updated 2020-06-29
# The script to do the ray tracing


import numpy as np
from sunRay import plasmaFreq as pfreq
from sunRay import densityModel as dm
from sunRay import scattering as scat 
#from sunRay import showPlot as SP   # matplotlib will slow down the run
from sunRay.parameters import c,c_r,R_S  # physics parameters
from sunRay.parameters import dev_u  # computation device
import sunRay.statisticalRays as raystat
import sunRay.record_properties as record_properties
import torch
import time
from tqdm import tqdm # for processing bar
import datetime

torch.set_default_tensor_type(torch.FloatTensor) # float is enough # float64 is overcare

def runRays(steps_N  = -1 , collect_N = 180, t_param = 20.0, photon_N = 10000,
            start_r = 1.75, start_theta = 0/180.0*np.pi,    start_phi  = 0/180.0*np.pi,
            f_ratio  = 1.1, ne_r = dm.parkerfit,    epsilon = 0.4, anis = 0.2,
            asym = 1.0, Te = 86.0, Scat_include = True, Show_param = True,
            Show_result_k = False, Show_result_r = False,  verb_out = False,
            sphere_gen = False, num_thread =4, early_cut= True ,dev_u = dev_u,
            save_npz = False, data_dir='./datatmp/',save_level=1,ignore_down=True,
            Absorb_include=True,dk_record=True):
    """
    name: runRays
    
    parameters:
        steps_N : number of the step # set as -1 to autoset
        collect_N : number of recorded step
        t_param : parameter of t step length 
            (larger t_parm corresponding to smaller dt)
        photon_N : number of photon
        start_r : start radius of the burst, in solar radii
        start_theta : in rad
        start_phi : in rad
        f_ratio : f/f_pe
        ne_r : density model used for this calculation 
        epsilon : fluctuation scale
        anis : the anisotropic parameter
        asym : asymetric scale
        Te : eV temperature in eV
        Scat_include : whether to consider the scattering  
        Show_param : Display the parameters
        Show_result_k : Show simulation result k
        verb_out : print message
        dev_u : device to use for the calculation
        save_npz [Bool] : whether to save the simulation result to file
        dir_npz : the directory for the npz data file 
        save_level : reduct the data to a certain level then save
        ignore_down: ignore the downward wave
        Louis: what is sphere_gen ?

    results:
        The t k and r of the ray-tracing result
    """
    print(photon_N)

    torch.set_num_threads(num_thread)
    # put variable in device
    start_r = torch.tensor([start_r])  
    PI = torch.acos(-torch.ones(1,device=dev_u))
    nu_e0 = 2.91e-6*ne_r(start_r)*20./Te**1.5     #Louis : I don't know what this physical parameter is
    nu_e = nu_e0
    photon_N_exist=photon_N
    
    # frequency of the wave
    freq0 = f_ratio * pfreq.omega_pe_r(ne_r,start_r.to(dev_u),dev_u=dev_u)/(2*PI)

    if verb_out:
        print('----------------------------------')
        print('Frequency : '+str(freq0.cpu().data.numpy()/1e6)[1:7]+'MHz')
        print('Compute with : '+str(dev_u))
        print('----------------------------------')

    #freq0 = torch.Tensor([freq0]).to(dev_u)

    # position of the photons
    rxx = start_r * torch.Tensor(np.sin(start_theta) * np.cos(start_phi) * np.ones(photon_N))
    ryy = start_r * torch.Tensor(np.sin(start_theta) * np.sin(start_phi) * np.ones(photon_N))
    rzz = start_r * torch.Tensor(np.cos(start_theta) * np.ones(photon_N))
    rr = start_r.to(dev_u) * torch.ones(photon_N).to(dev_u)
    rr_cur = rr # rr_cur [current rr for for loop]
    r_vec = torch.stack((rxx,ryy,rzz),0).to(dev_u)

    omega0 = freq0*(2*PI)
    #nu_s0 = scat.nuScattering(rr,omega0,epsilon,ne_r,dev_u=dev_u)
    nu_s0 = scat.nuScatterKrupar(rr,omega0,epsilon,ne_r,dev_u=dev_u)

    #if Show_param:
    #    SP.showParameters(ne_r,omega0,epsilon)  

    # wave-vector of the photons
    kc0 = torch.sqrt(omega0**2. - pfreq.omega_pe_r(ne_r,rr,dev_u=dev_u)**2.)

    t_dist = np.zeros(photon_N)   #a variable defined by Louis
    dict_dist = {'time': torch.zeros(photon_N), 'x': torch.zeros(photon_N), 'y': torch.zeros(photon_N), 'z': torch.zeros(photon_N)}
    
    if sphere_gen:
        k_theta = torch.Tensor(np.random.uniform(low=-np.pi/2 + 1e-4 ,
                                high= np.pi/2 ,size=photon_N)).to(dev_u) # k_z > 0 
        k_mu0   = torch.cos(k_theta)
        k_phi0  = torch.Tensor(np.random.uniform(low=0 ,
                                high= 2*np.pi, size=photon_N)).to(dev_u) # phi in all dir

        kxx_k = kc0 * torch.sqrt(1-k_mu0**2.) * torch.cos(k_phi0)
        kyy_k = kc0 * torch.sqrt(1-k_mu0**2.) * torch.sin(k_phi0)
        kzz_k = kc0 * k_mu0
        k_vec = torch.stack((kxx_k,kyy_k,kzz_k),0).to(dev_u)
    else:
        # generate in xyz
        k_vec_tmp = torch.randn(3,photon_N).to(dev_u)
        k_vec = kc0 * k_vec_tmp/torch.sqrt(torch.sum(k_vec_tmp.pow(2),axis=0))
        # ignore downward (r k not same direction)
        idx_select = torch.nonzero(torch.sum(r_vec*k_vec,axis=0)<0,as_tuple=False)
        if ignore_down:
            k_vec[:,idx_select] = -k_vec[:,idx_select] 

    r_vec_start = r_vec
    k_vec_start = k_vec

    
    kc = torch.sqrt(torch.sum(k_vec.pow(2),axis=0))
    kc_cur = kc
    

    # Detach from the previous compute graph
    # before record steps for diff
    domega_pe_dxyz = pfreq.domega_dxyz_1d(ne_r,r_vec.detach(),dev_u=dev_u)

    Exp_size = 1.25*30./(freq0/1e6)
    dt0 = 0.01*Exp_size/c_r
    tau = torch.zeros(rr_cur.shape).to(dev_u)
    dk_inte_refr = torch.zeros(rr_cur.shape).to(dev_u)
    dk_inte_scat = torch.zeros(rr_cur.shape).to(dev_u)
    if dk_record:
        dkx_inte_refr = torch.zeros(rr_cur.shape).to(dev_u)
        dkx_inte_scat = torch.zeros(rr_cur.shape).to(dev_u)
        dky_inte_refr = torch.zeros(rr_cur.shape).to(dev_u)
        dky_inte_scat = torch.zeros(rr_cur.shape).to(dev_u)
        dkz_inte_refr = torch.zeros(rr_cur.shape).to(dev_u)
        dkz_inte_scat = torch.zeros(rr_cur.shape).to(dev_u)
    
    # a function to find the 1/1e4 (should be 1/e3) small element in the array
    find_small_1e3 = lambda arr:  torch.sort(arr)[0][int(photon_N*1e-3)]

    if steps_N == -1:
        dt_dr0  = find_small_1e3(rr_cur/omega0*kc_cur)/t_param
        dt_nu0  = find_small_1e3(1/(nu_s0)) 
        dt_nue0  = 1/nu_e0
        steps_N = (1.5*4.605/nu_e0/(1.-1./f_ratio**2)**0.5 + 
            15/c_r)*(1/dt_dr0+1.5/dt_nu0+1/dt_nue0+25)*(0.3+(anis*4)) + 8192  #(0.1+(anis**0.5)) 
        #Louis: what is this calculation for steps_N based on???
        #the larger start_r, the larger steps_N 
        #steps_N = torch.tensor([50000])  
        if verb_out:
            print("Refraction dt : "+str(1/dt_dr0.cpu().numpy()))
            print("Scattering dt : "+str(1/dt_nu0.cpu().numpy()))
            print("Absorb Col    : "+str(1/dt_nue0.cpu().numpy()[0]))
            print("Absorb  t     : "+str((1.5*4.605/nu_e0/(1.-1./f_ratio**2)**0.5)[0].cpu().numpy()))
    else :
        dt_dr0  = find_small_1e3(rr_cur/omega0*kc_cur)/t_param
        dt_nu0  = find_small_1e3(1/(nu_s0)) 
        dt_nue0  = 1/nu_e0
        steps_N = torch.tensor([steps_N])       
    
    # collect the variables of the simulation
    # collect to CPU (GPU mem is expensive)
    collectPoints = np.round(np.linspace(0,steps_N-1,collect_N))
    #collectPoints = np.arange(0, collect_N, 1)  #Louis : I added this so that we can record more steps




    r_vec_collect = torch.zeros(collect_N,3,photon_N).to(torch.device('cpu'))-1   #orginally has 'cpu' as argument 
    k_vec_collect = torch.zeros(collect_N,3,photon_N).to(torch.device('cpu'))-1   #orginally has 'cpu' as argument 
    tau_collect = torch.zeros(collect_N,photon_N).to(torch.device('cpu'))-1       #orginally has 'cpu' as argument 
    t_collect = torch.zeros(collect_N).to(torch.device('cpu'))-1                  #orginally has 'cpu' as argument 
    dk_refr_collect = torch.zeros(collect_N,photon_N).to(torch.device('cpu'))-1   #orginally has 'cpu' as argument 
    dk_scat_collect = torch.zeros(collect_N,photon_N).to(torch.device('cpu'))-1   #orginally has 'cpu' as argument 

    
    if dk_record:
        dkx_refr_collect = torch.zeros(collect_N,photon_N).to(dev_u)
        dkx_scat_collect = torch.zeros(collect_N,photon_N).to(dev_u)
        dky_refr_collect = torch.zeros(collect_N,photon_N).to(dev_u)
        dky_scat_collect = torch.zeros(collect_N,photon_N).to(dev_u)
        dkz_refr_collect = torch.zeros(collect_N,photon_N).to(dev_u)
        dkz_scat_collect = torch.zeros(collect_N,photon_N).to(dev_u)

    else:
        (dkx_refr_collect,dky_refr_collect,dkz_refr_collect,
        dkx_scat_collect,dkx_scat_collect,dkx_scat_collect)=(0,0,0,0,0,0)
        
    idx_collect  =  0
    t_current = 0

    
    time.sleep(0.5)  #this causes problems for some reason
    
    pbar=tqdm(np.arange(steps_N))




    DT = []
    T = []
    print('The fixed timestep is:',(40/c_r/(steps_N/collect_N))[0])
    # the big loop
    for idx_step in (pbar if verb_out else np.arange(steps_N)): #show process bar
    #for idx_step in (np.arange(steps_N)):
            
        # dispersion relation reform
        omega = torch.sqrt(pfreq.omega_pe_r(ne_r,rr_cur,dev_u=dev_u)**2 + kc_cur**2)
        freq_pe = omega/(2*PI)

        nu_s = scat.nuScattering(rr_cur,omega,epsilon,ne_r,dev_u=dev_u)
        nu_s = nu_s*(nu_s<nu_s0)+nu_s0*(~(nu_s<nu_s0)) # use the smaller nu_s

        # detach is essential for remove the previous calc map
        domega_pe_dxyz = pfreq.domega_dxyz_1d(ne_r,r_vec.detach(),dev_u=dev_u)
        domega_pe_dr = torch.sqrt(torch.sum(domega_pe_dxyz.pow(2),axis=0))

        with torch.no_grad(): # no autograd in following calc
            dr_vec = c_r/omega.repeat(3,1) * k_vec    

            # component of r and k vector at current step
            rr_cur = torch.sqrt(torch.sum(r_vec.pow(2),axis=0))
            kc_cur = torch.sqrt(torch.sum(k_vec.pow(2),axis=0))
            rx_cur,ry_cur,rz_cur = r_vec[0,:],r_vec[1,:],r_vec[2,:]
            kx_cur,ky_cur,kz_cur = k_vec[0,:],k_vec[1,:],k_vec[2,:]


            #code from Louis
            #============================================================================================================================================
            #if (idx_step%10000==0):
                #print('Should we stop?:', find_small_1e3(rr_cur)>111)
            dict_dist = record_properties.record_properties2(215, dict_dist, rr_cur, rx_cur, ry_cur, rz_cur, t_current)
            #t_dist = record_properties.record_properties_thetarng(111, t_dist, 0, 10, rr_cur, rx_cur, ry_cur, rz_cur, t_current)

            #============================================================================================================================================


            # dynamic time step
            # for large particle size, use a part to estimate
            if photon_N>10001: 
                dt_ref = find_small_1e3((torch.abs(kc_cur/ (domega_pe_dr*c_r)/t_param))[0:10000]) # t step
                dt_dr  = find_small_1e3((rr_cur/omega0*kc_cur)[0:10000])/t_param
                dt_nu  = find_small_1e3((0.1/(nu_s))[0:10000]) 
            else:    
                dt_ref = find_small_1e3(torch.abs(kc_cur/ (domega_pe_dr*c_r)/t_param)) # t step
                dt_dr  = find_small_1e3(rr_cur/omega0*kc_cur)/t_param
                dt_nu  = find_small_1e3(0.1/(nu_s)) 
            dt_fix=(40/c_r/(steps_N/collect_N))[0]   #Louis: this will be the same everytime, so if you want 
                                                     #constant timestepping you want dt to be this at every
                                                     #timestep
            
            # make sure most of the photons have proper dt 
            dt = torch.Tensor([np.nanmin([dt_nu,dt_ref,dt_dr,dt0,dt_fix])]).to(dev_u)
            dt_idx = np.argmin(([dt_nu,dt_ref,dt_dr,dt0,dt_fix]))
            #dt = torch.tensor([0.01])
            #print(dt_nu,dt_ref,dt_dr,dt0,dt_fix)
            #dt = torch.Tensor([np.nanmin([0.06])]).to(dev_u)   #Louis

            DT.append([dt, dt_idx])
            T.append(t_current)


            #if verb_out :
            #    pbar.set_postfix({'dt': [dt.cpu().numpy()[0],
            #                             dt_ref.cpu(),
            #                             dt_dr.cpu(),
            #                             dt_nu.cpu(),
            #                             dt0.cpu()]})

            g0 = torch.sqrt(nu_s*kc_cur**2)

            # random vec for wave scattering  # [3*N] normal distribution
            W_vec = -torch.randn(r_vec.shape,device=dev_u) * torch.sqrt(dt) 
            #W_vec = torch.randn(r_vec.shape).to(dev_u) * torch.sqrt(dt)   # slow
            Wx,Wy,Wz = W_vec[2,:],W_vec[1,:],W_vec[0,:]

            # photon position in spherical coordinates
            # (rx,ry,rz) is the direction of anisotropic tubulence
            fi = torch.atan2(ry_cur,rx_cur)
            costheta = rz_cur/rr_cur
            sintheta = torch.sqrt(1-costheta**2)
            if Scat_include:

                # rotate the k vec into the r-z coordinate
                kcx = - kx_cur*torch.sin(fi) + ky_cur*torch.cos(fi) 
                kcy = (- kx_cur*costheta*torch.cos(fi) 
                    - ky_cur*costheta*torch.sin(fi) + kz_cur*sintheta) 
                kcz = (  kx_cur*sintheta*torch.cos(fi) 
                    + ky_cur*sintheta*torch.sin(fi) + kz_cur*costheta)

                kw     =  Wx*kcx+Wy*kcy+Wz*kcz*anis
                Akc    = torch.sqrt(kcx*kcx+kcy*kcy+kcz*kcz*(anis**2))
                z_asym = (asym*(kcz > 0.0) + (2.0-asym)*(~(kcz>0.))) * (kc_cur/Akc)**2

                A_perp = (nu_s*z_asym* kc_cur /(Akc**3) *
                    (-(1+anis**2)*Akc**2+3*anis**2 *(anis**2-1)*kcz**2) *anis)
                A_par  = (nu_s*z_asym* kc_cur /(Akc**3) *
                    ((-3*anis**4+anis**2)*(Akc**2)+3*anis**4 * (anis**2-1)*kcz**2)*anis)
                A_g0   = g0*torch.sqrt(z_asym*anis)

                dkx_tmp = A_perp*kcx*dt + A_g0*(Wx-kcx*kw/Akc**2)
                dky_tmp = A_perp*kcy*dt + A_g0*(Wy-kcy*kw/Akc**2)
                dkz_tmp = A_par *kcz*dt + A_g0*(Wz-kcz*kw*anis/Akc**2)*anis

                kcx = kcx + dkx_tmp
                kcy = kcy + dky_tmp
                kcz = kcz + dkz_tmp

                dk_inte_scat = dk_inte_scat + torch.sqrt((dkx_tmp/omega0)**2+
                            (dky_tmp/omega0)**2+(dkz_tmp/omega0)**2)*c_r
                
                # rotate back to normal coordinate
                kx_cur_new = (-kcx*torch.sin(fi) 
                    -kcy*costheta*torch.cos(fi) +kcz*sintheta*torch.cos(fi) )
                ky_cur_new = ( kcx*torch.cos(fi) 
                    -kcy*costheta*torch.sin(fi) +kcz*sintheta*torch.sin(fi) )
                kz_cur_new =  kcy*sintheta+kcz*costheta

                
                if dk_record:
                    dkx_inte_scat = dkx_inte_scat + ((kx_cur_new-kx_cur)/omega0)*c_r
                    dky_inte_scat = dky_inte_scat + ((ky_cur_new-ky_cur)/omega0)*c_r
                    dkz_inte_scat = dkz_inte_scat + ((kz_cur_new-kz_cur)/omega0)*c_r
                
                kx_cur,ky_cur,kz_cur= kx_cur_new,ky_cur_new,kz_cur_new

            r_vec = torch.stack((rx_cur,ry_cur,rz_cur),0)
            k_vec = torch.stack((kx_cur,ky_cur,kz_cur),0)

            rr_cur = torch.sqrt(torch.sum(r_vec.pow(2),axis=0))
            kc_cur = torch.sqrt(torch.sum(k_vec.pow(2),axis=0))

            # re-normalize  # to keep |k_vec| stable
            kc_norm = torch.sqrt(kx_cur**2 + ky_cur**2 + kz_cur**2)
            k_vec = k_vec * kc_cur.repeat(3,1)/ kc_norm.repeat(3,1)

            # k step forward  # refraction
            dk_xyz_dt = ((pfreq.omega_pe_r(ne_r,rr_cur,dev_u=dev_u)/omega).repeat(3,1)   
                        * domega_pe_dxyz) * c_r
            k_vec = k_vec - dk_xyz_dt * dt

            dk_inte_refr = dk_inte_refr + torch.sqrt(torch.sum(
                    (dk_xyz_dt*dt/omega0).pow(2),axis=0))*c_r

            if dk_record:
                dkx_inte_refr = dkx_inte_refr + (dk_xyz_dt*dt/omega0)[0,:]*c_r
                dky_inte_refr = dky_inte_refr + (dk_xyz_dt*dt/omega0)[1,:]*c_r
                dkz_inte_refr = dkz_inte_refr + (dk_xyz_dt*dt/omega0)[2,:]*c_r
            
            # r step forward
            r_vec = r_vec + dr_vec * dt

            # update abs after vec change        
            rr_cur = torch.sqrt(torch.sum(r_vec.pow(2),axis=0))
            kc_cur = torch.sqrt(torch.sum(k_vec.pow(2),axis=0))

            # re-normalize  # to keep omega stable
            kc_refresh = (torch.sqrt(omega**2-pfreq.omega_pe_r(ne_r,rr_cur,dev_u=dev_u)**2)
                /torch.sqrt(torch.sum(k_vec.pow(2),axis=0)))
            k_vec = k_vec * kc_refresh.repeat(3,1)

            nu_e = (2.91e-6 * ne_r(rr_cur) * 20. / Te**1.5
                *pfreq.omega_pe_r(ne_r, rr_cur,dev_u=dev_u)**2/omega**2)
            
            tau = tau + nu_e*dt

            rr_cur = torch.sqrt(torch.sum(r_vec.pow(2),axis=0))
            kc_cur = torch.sqrt(torch.sum(k_vec.pow(2),axis=0))

            # absorb the photon with large optical depth(set as NaN)
            # 11.513 -> I=1e-5 for very low frequency ratio
            # 9.210  -> I=1e-4
            # 6.908  -> I=1e-3
            # 4.605  -> I=1e-2
            # remove every 128 steps
            if (idx_step%128==0) and Absorb_include:
                idx_absorb = torch.nonzero(tau>11.513,as_tuple=False)
                r_vec[:,idx_absorb] = r_vec[:,idx_absorb]*torch.Tensor([np.nan]).to(dev_u) 
                k_vec[:,idx_absorb] = k_vec[:,idx_absorb]*torch.Tensor([np.nan]).to(dev_u) 
                rr_cur[idx_absorb] =  rr_cur[idx_absorb]*torch.Tensor([np.nan]).to(dev_u) 
                kc_cur[idx_absorb] =  kc_cur[idx_absorb]*torch.Tensor([np.nan]).to(dev_u)
                tau[idx_absorb] = tau[idx_absorb]*torch.Tensor([np.nan]).to(dev_u)

                # remove [tail and back propagation]
                absorb_tail=False
                if absorb_tail:
                    idx_absorb2 = torch.nonzero( ((torch.sum(r_vec*k_vec,axis=0)/(rr_cur*kc_cur))<0.01) & 
                                                (rr_cur < find_small_1e3(rr_cur)),
                                                as_tuple=False)
                
                    r_vec[:,idx_absorb2] = r_vec[:,idx_absorb2]*torch.Tensor([np.nan]).to(dev_u) 
                    k_vec[:,idx_absorb2] = k_vec[:,idx_absorb2]*torch.Tensor([np.nan]).to(dev_u) 
                    rr_cur[idx_absorb2] =  rr_cur[idx_absorb2]*torch.Tensor([np.nan]).to(dev_u) 
                    kc_cur[idx_absorb2] =  kc_cur[idx_absorb2]*torch.Tensor([np.nan]).to(dev_u)
                    tau[idx_absorb2] = tau[idx_absorb2]*torch.Tensor([np.nan]).to(dev_u)
                    
                # change size every 1024 step [drop removed photons]
                change_size=False
                if (idx_step%2048==0):
                    if change_size:
                        idx_exist = ~torch.isnan(tau)
                        tau=tau[idx_exist]
                        tau_collect = tau_collect[:,idx_exist]
                        k_vec=k_vec[:,idx_exist]
                        r_vec=r_vec[:,idx_exist]
                        k_vec_collect=k_vec_collect[:,:,idx_exist]
                        r_vec_collect=r_vec_collect[:,:,idx_exist]
                        dk_refr_collect=dk_refr_collect[:,idx_exist]
                        dk_scat_collect=dk_scat_collect[:,idx_exist]

                        if dk_record:
                            dkx_refr_collect=dkx_refr_collect[:,idx_exist]
                            dky_refr_collect=dky_refr_collect[:,idx_exist]
                            dkz_refr_collect=dkz_refr_collect[:,idx_exist]
                            dkx_scat_collect=dkx_scat_collect[:,idx_exist]
                            dky_scat_collect=dky_scat_collect[:,idx_exist]
                            dkz_scat_collect=dkz_scat_collect[:,idx_exist]
                    else:   
                        photon_N_exist=tau[~torch.isnan(tau)].shape[0]
                        find_small_1e3 = lambda arr:  torch.sort(arr)[0][int(photon_N_exist*1e-3)]


        t_current = t_current + dt
        #temp = True #Louis
        if idx_step in collectPoints:
        #if idx_step % 50 == 0:
        #if temp:   #Louis
            #print(idx_collect)
            print('Collect at:', idx_step)
            t_collect[idx_collect] = t_current
            r_vec_collect[idx_collect,:,:] = r_vec.cpu()
            k_vec_collect[idx_collect,:,:] = k_vec.cpu()
            tau_collect[idx_collect,:] = tau.cpu()
            dk_refr_collect[idx_collect,:] = dk_inte_refr.cpu()
            dk_scat_collect[idx_collect,:] = dk_inte_scat.cpu()
            if dk_record:
                dkx_refr_collect[idx_collect,:] = dkx_inte_refr.cpu()
                dky_refr_collect[idx_collect,:] = dky_inte_refr.cpu()
                dkz_refr_collect[idx_collect,:] = dkz_inte_refr.cpu()
                dkx_scat_collect[idx_collect,:] = dkx_inte_scat.cpu()
                dky_scat_collect[idx_collect,:] = dky_inte_scat.cpu()
                dkz_scat_collect[idx_collect,:] = dkz_inte_scat.cpu()
                ## to be finished
            if idx_collect < collect_N -1:  #Louis: I've added this condition!
                idx_collect = idx_collect +1
            if verb_out==2: # print out the process
                print('F_pe:'+'{:.3f}'.format(np.mean(
                    (pfreq.omega_pe_r(ne_r,torch.mean(rr_cur),dev_u=dev_u)/2/PI/1e6).cpu().data.numpy()))+
                    ' |  R:'+'{:.3f}'.format(torch.mean(rr_cur).cpu().data.numpy())+
                    ' |  Ne_r:'+'{:.3f}'.format(ne_r(torch.mean(rr_cur)).cpu().data.numpy())+
                    ' |  nu_s: ' +  '{:.3f}'.format(torch.mean(0.1/nu_s).cpu().data.numpy())+
                    ' |  F_ratio: ' +  '{:.3f}'.format(torch.mean(omega0/omega).cpu().data.numpy()))
            
        if early_cut and (idx_step>5000):
            if find_small_1e3(rr_cur)>220: # all out of 1AU, so the simulation stops once all the photons are further away than 1 AU
                #print('Has there something been removed?:', r_vec_collect.shape)
                # remove nans
                idx_exist = ~torch.isnan(tau)
                tau=tau[idx_exist]
                tau_collect = tau_collect[:,idx_exist]
                k_vec=k_vec[:,idx_exist]
                r_vec=r_vec[:,idx_exist]
                k_vec_collect=k_vec_collect[:,:,idx_exist]
                r_vec_collect=r_vec_collect[:,:,idx_exist]
                dk_refr_collect=dk_refr_collect[:,idx_exist]
                dk_scat_collect=dk_scat_collect[:,idx_exist]

                if dk_record:
                    dkx_refr_collect=dkx_refr_collect[:,idx_exist]
                    dky_refr_collect=dky_refr_collect[:,idx_exist]
                    dkz_refr_collect=dkz_refr_collect[:,idx_exist]
                    dkx_scat_collect=dkx_scat_collect[:,idx_exist]
                    dky_scat_collect=dky_scat_collect[:,idx_exist]
                    dkz_scat_collect=dkz_scat_collect[:,idx_exist]
                                    
                photon_N_exist=tau.shape[0]
                find_small_1e3 = lambda arr:  torch.sort(arr)[0][int(photon_N_exist*1e-3)]
                
                
                final_collect = idx_collect
                # cut
                t_collect = t_collect[0:final_collect]
                
                r_vec_collect = r_vec_collect[0:final_collect,:,:] 
                k_vec_collect = k_vec_collect[0:final_collect,:,:] 
                tau_collect   = tau_collect[0:final_collect,:] 
                dk_refr_collect =dk_refr_collect[0:final_collect,:] 
                dk_scat_collect =dk_scat_collect[0:final_collect,:] 
                
                if dk_record:
                    dkx_refr_collect =dkx_refr_collect[0:final_collect,:] 
                    dkx_scat_collect =dkx_scat_collect[0:final_collect,:] 
                    dky_refr_collect =dky_refr_collect[0:final_collect,:] 
                    dky_scat_collect =dky_scat_collect[0:final_collect,:] 
                    dkz_refr_collect =dkz_refr_collect[0:final_collect,:] 
                    dkz_scat_collect =dkz_scat_collect[0:final_collect,:] 

                #print('We stopped at step:', idx_step)
                #temp = np.copy(rr_cur)
                #print(rr_cur, temp)
                collect_N = final_collect

                break # stop the loop
            


    t_collect_local = t_collect.cpu().data.numpy()
    r_vec_collect_local  = r_vec_collect.cpu().data.numpy()
    k_vec_collect_local  = k_vec_collect.cpu().data.numpy()
    tau_collect_local = tau_collect.cpu().data.numpy()

    #plt.figure(1)
    #plt.plot(t_collect_local, r_vec_collect_local[:,0,0])
    #plt.figure(2)
    #plt.plot( r_vec_collect_local[:,0,0], r_vec_collect_local[:,1,0])
    print('Traced final t : '+str(t_collect_local[-1])+' s')

    #if Show_result_r:
    #    SP.showResultR(r_vec_collect_local)

    if save_npz:
         # save the data to npz file
        if save_level == 0:
            import datetime
            t_stamp=str(datetime.datetime.now()).replace(' ','_').replace('-','').replace(':','')[0:15]
            np.savez_compressed(data_dir+'RUN_[eps'+str(np.round(epsilon,5)) +
                ']_[alpha'+str(np.round(anis,5))+']_R'+str(np.round(f_ratio,5))+'_'+t_stamp+'.lv1.npz', 
                steps_N  = steps_N, 
                collect_N = collect_N, photon_N = photon_N, start_r = start_r, 
                start_theta = start_theta, start_phi  = start_phi, 
                f_ratio  = f_ratio, epsilon = epsilon , anis = anis, asym = asym,
                omega0=omega0.cpu(), freq0=freq0.cpu(),
                t_collect=t_collect.cpu(), tau=tau.cpu(),
                r_vec_collect_local=r_vec_collect_local,
                k_vec_collect_local=k_vec_collect_local,
                tau_collect_local = tau_collect_local)
        if save_level == 1:

            (r_vec_stat_avail,k_vec_stat_avail,t_reach_stat_avail,tau_stat_avail,
                r_vec_0, k_vec_0) =  raystat.reduct_lv1(
                    photon_N,r_vec_collect_local,k_vec_collect_local,
                    t_collect,tau_collect_local,omega0,num_t_bins=60)
            
            if dk_record:
                # for dk of refraction
                (r_vec_stat_avail,k_vec_stat_avail,t_reach_stat_avail,dk_refr_avail,
                                r_vec_0, k_vec_0) =  raystat.reduct_lv1(
                                    photon_N,r_vec_collect_local,k_vec_collect_local,
                                    t_collect,dk_refr_collect,omega0,num_t_bins=60)
                # for dk of scattering
                (r_vec_stat_avail,k_vec_stat_avail,t_reach_stat_avail,dk_scat_avail,
                                r_vec_0, k_vec_0) =  raystat.reduct_lv1(
                                    photon_N,r_vec_collect_local,k_vec_collect_local,
                                    t_collect,dk_scat_collect,omega0,num_t_bins=60)
                # for dk_vec scat
                (r_vec_stat_avail,k_vec_stat_avail,t_reach_stat_avail,dkx_scat_avail,
                                r_vec_0, k_vec_0) =  raystat.reduct_lv1(
                                    photon_N,r_vec_collect_local,k_vec_collect_local,
                                    t_collect,dkx_scat_collect,omega0,num_t_bins=60)
                (r_vec_stat_avail,k_vec_stat_avail,t_reach_stat_avail,dky_scat_avail,
                                r_vec_0, k_vec_0) =  raystat.reduct_lv1(
                                    photon_N,r_vec_collect_local,k_vec_collect_local,
                                    t_collect,dky_scat_collect,omega0,num_t_bins=60)
                (r_vec_stat_avail,k_vec_stat_avail,t_reach_stat_avail,dkz_scat_avail,
                                r_vec_0, k_vec_0) =  raystat.reduct_lv1(
                                    photon_N,r_vec_collect_local,k_vec_collect_local,
                                    t_collect,dkz_scat_collect,omega0,num_t_bins=60)
                
                # for dk_vec refr
                (r_vec_stat_avail,k_vec_stat_avail,t_reach_stat_avail,dkx_refr_avail,
                                r_vec_0, k_vec_0) =  raystat.reduct_lv1(
                                    photon_N,r_vec_collect_local,k_vec_collect_local,
                                    t_collect,dkx_refr_collect,omega0,num_t_bins=60)
                (r_vec_stat_avail,k_vec_stat_avail,t_reach_stat_avail,dky_refr_avail,
                                r_vec_0, k_vec_0) =  raystat.reduct_lv1(
                                    photon_N,r_vec_collect_local,k_vec_collect_local,
                                    t_collect,dky_refr_collect,omega0,num_t_bins=60)
                (r_vec_stat_avail,k_vec_stat_avail,t_reach_stat_avail,dkz_refr_avail,
                                r_vec_0, k_vec_0) =  raystat.reduct_lv1(
                                    photon_N,r_vec_collect_local,k_vec_collect_local,
                                    t_collect,dkz_refr_collect,omega0,num_t_bins=60)
            else:
                dk_refr_avail,dk_scat_avail = 0,0
                dkx_refr_avail,dky_refr_avail,dkz_refr_avail,dkx_scat_avail,dky_scat_avail,dkz_scat_avail=0,0,0,0,0,0
            
            import datetime
            t_stamp=str(datetime.datetime.now()).replace(' ','_').replace('-','').replace(':','')[0:15]
            np.savez_compressed(data_dir+'RUN_[eps'+str(np.round(epsilon,5)) +
                ']_[alpha'+str(np.round(anis,5))+']_R'+str(np.round(f_ratio,5))+'_'+t_stamp+'.lv1.npz', 
                steps_N  = steps_N, 
                collect_N = collect_N, photon_N = photon_N, start_r = start_r, 
                start_theta = start_theta, start_phi  = start_phi, 
                f_ratio  = f_ratio, epsilon = epsilon , anis = anis, asym = asym,
                omega0=omega0.cpu(), freq0=freq0.cpu(),
                r_vec_stat_avail=r_vec_stat_avail,k_vec_stat_avail=k_vec_stat_avail,
                t_reach_stat_avail=t_reach_stat_avail,
                tau_stat_avail=tau_stat_avail,r_vec_0=r_vec_0, k_vec_0=k_vec_0,
                dk_refr_avail=dk_refr_avail,dk_scat_avail=dk_scat_avail,dkx_refr_avail=dkx_refr_avail,
                dky_refr_avail=dky_refr_avail,
                dkz_refr_avail=dkz_refr_avail,
                dkx_scat_avail=dkx_scat_avail,
                dky_scat_avail=dky_scat_avail,
                dkz_scat_avail=dkz_scat_avail)

    dict_dist['time'], dict_dist['x'], dict_dist['y'], dict_dist['z'] = dict_dist['time'].to('cpu'), dict_dist['x'].to('cpu'), dict_dist['y'].to('cpu'), dict_dist['z'].to('cpu') 
            
    return (steps_N  ,  collect_N,  photon_N, start_r,  start_theta, start_phi,  f_ratio, 
            epsilon ,  anis, asym,  omega0.cpu(), freq0.cpu(), t_collect.cpu(), tau.cpu(),
            r_vec_collect_local,  k_vec_collect_local,  tau_collect_local,
            dk_refr_collect, dk_scat_collect, dict_dist, DT, T)

