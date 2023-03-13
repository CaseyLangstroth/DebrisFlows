# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 13:06:40 2022

@author: Scott
"""

       #  # if hits the end of the domain append the full distance traveled to 
       #  # cummulative celluslar distance 
       #  if tdst<dist_downstream:
       #      dst_down=np.append(dst_down,dist_downstream)
        
       #  # build the profile of the wedge e.g. how thick is the sediment
       #  sed_thick=np.tan(beta)*dst_down
        
       #  dl=len(sed_profile)
       # # Ased=np.zeros(dl)
       #  Ased=np.zeros(dl)
       #  print('change to river width')
       #  VBi=VBw_ras[xi,yi]
        
       #  for k in np.arange(0,dl-1):  
       #    sed_thick[k+1]=np.trapz([sed_profile[k],sed_profile[k+1]],[dst_down[k],dst_down[k+1]])
            
        # old version to compute area and convert to sediment volume
        #     Ased[k+1]=np.trapz([sed_profile[k],sed_profile[k+1]],[dst_down[k],dst_down[k+1]])

        #     #print(k)
        #     #print(sed_profile[k])
        # vol_sed_down=Ased*VBi
        # #Input to River*\
        # DownValley_Input[xi,yi] = DownValley_Input[xi,yi]+vol_sed_down
        #DVi=DownValley_Input.copy()
    #Total Volume Input To River