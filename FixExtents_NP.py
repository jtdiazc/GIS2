import numpy as np
import os

np_dir=r'\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\Accretion\20211223'

#Import template array

temp_np_nam="SLR_1_1"
temp_np=np.load(os.path.join(np_dir,temp_np_nam+".npy"))

#Import array to be fixed
tbf_np_nam="MTL"
tbf_np=np.load(os.path.join(np_dir,tbf_np_nam+".npy"))

tbf_np[np.isnan(temp_np)]=np.nan

np.save(os.path.join(np_dir,tbf_np_nam+'_fix.npy'),tbf_np)
