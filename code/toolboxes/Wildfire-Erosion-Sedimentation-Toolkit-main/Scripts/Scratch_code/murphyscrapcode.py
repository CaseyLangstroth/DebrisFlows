# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 14:53:14 2023

@author: Scott
"""

        '''
        # Compute initial thickness and volumes 
        # calculate thickness along wedge assuming uniform width
        sed_profile=np.tan(beta)*dst_up
        sed_vol=np.zeros(len(sed_profile)-1)
        dl=len(sed_vol)
        
        for k in np.arange(0,dl):
            sed_vol[k]=np.trapz([sed_profile[k+1],sed_profile[k]],x=[dst_up[k+1],dst_up[k]])
        
        # conservation of mass checking
        tvol_init=np.max(np.cumsum(sed_vol)*b)
            # /dist_step(k)])
        
        
        # # get distance between each cell
        # dx_step1=np.abs(np.diff(dst_up))
        # print('change 10 here to raster resolution')
        # dx_step=np.insert(dx_step1,0,10.0)
        
        # get channel widths in area of interest
        chanwidth=RW_ras[xi,yi]
        
        # if overdomain:
        #     chanwidth=np.append(chanwidth,chanwidth[-1])
        # # multiply them together to get an inital total volume per cell
        # vol_init=b*dx_step*sed_thick
        
        # tvol_init=np.max(np.cumsum(vol_init))
        
        
        # for k in np.arange(0,dl-1):
        #     sed_thick[k+1]=np.trapz([sed_profile[k],sed_profile[k+1]],[dst_down[k],dst_down[k+1]])/dist_step(k)
        
        
        # get channel widths in area of interest
        #chanwidth=RW_ras[xi,yi]
        # # get distance between each cell this accounts for diagonal steps
        # dstep=np.diff(dst_down)
        
        # # append zero to front of array for the inital cell
        # dstep2=np.insert(dstep,0,0)
        
        # reconstuct the total distance with starting value being the largest
        

        
        

            
            
        # sed_thick=np.tan(beta)*np.flipud(dst_down)
        # dist_step=np.diff(dst_down)
        # dl=len(dist_step)