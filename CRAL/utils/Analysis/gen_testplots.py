import yt_funcs
import cfuncs

Dir_sim='/ptmp/mpa/sraghu/noneq/' 

print("Dir_sim :",Dir_sim)

simdir_1='/u/sraghu/sim/SHUB/b256/sn_eq'
simdir_2='/u/sraghu/sim/SHUB/b256/sn_neq'
simdir_3='/u/sraghu/sim/SHUB/b256/sn_neq_18'
simdir_4='/u/sraghu/sim/SHUB/b256/sn_neq_lightstar'

simtyp_arr    = ['sn_eq', 'sn_neq', 'sn_neq_highres', 'sn_neq_lightstar']
#prop_arr      = ['sfr', 'hmf' ,  'dens_proj', 'temp_proj' , 'phaseplot']#,  'globalsink_vs_z',  'haloprops']

prop_arr = ['temp_proj']

simdir_arr     = [simdir_1, simdir_2, simdir_3, simdir_4]#, simdir_4]
cbar_arr       = ['m'     ,'r'      , 'g'     ,  'b'    ]#, 'b'     , 'y'     ]#, 'm'     ]#, 'm', 'c', 'k']
linestyle_arr  = ['-'     , '-'     , '-'     ,  '-']#, '-'     , '-'     ]
mbar_arr       = ["o"     , "^"     , "s"     ,  'h']#, "h"     , "D"     ]
sfr_outnum_arr = [11      , 11      , 11      ,  '11']#,8        ,9        ]#,8        ]
             #z8  #z9
outnum_1 = [ 10  , 9 ]
outnum_2 = [ 11  , 9 ]
outnum_3 = [ 11  , 9 ]
outnum_4 = [ 10  , 9 ]
Nout     = len(outnum_1)

outnum_arr = [ outnum_1, outnum_2, outnum_3, outnum_4 ]


for prop in prop_arr :
  print('######################## prop :', prop, "##########################")
  if(prop=='sfr') :
    #compare SFR
    plotanal=False
    plotdata=True
    overwrite=False
    #write star files
    iout=0
    for simdir in simdir_arr:
      outnum = sfr_outnum_arr[iout]
      yt_funcs.write_starfile(simdir, outnum, overwrite)
      iout = iout+1
    # plot sfr  
    yt_funcs.compare_sfr(simdir_arr, simtyp_arr ,sfr_outnum_arr, cbar_arr, linestyle_arr ,plotanal, plotdata )
  if(prop=='hmf') :
    # write halo catalog 
    overwrite=False
    for iout in range(0,Nout) :
      isim=0
      hmf_outnum_arr = []
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]  
        hmf_outnum_arr.append(outnum)
        yt_funcs.create_halo_catalog(simdir, outnum, 0, overwrite)
        isim = isim + 1
      fitfunc_arr = []
      print('hmf_outnum_arr : ', hmf_outnum_arr)
      yt_funcs.compare_hmf(simdir_arr, simtyp_arr, cbar_arr, mbar_arr, hmf_outnum_arr, fitfunc_arr)
  if(prop=='dens_proj') :
    boxlen_comov=6.77
    overwrite=False
    for iout in range(0,Nout) :
      isim=0
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        yt_funcs.plot_proj_gasdens(simdir, outnum, simtyp, boxlen_comov, overwrite)
        isim = isim + 1
  if(prop=='temp_proj') :
    boxlen_comov=6.77
    overwrite=False
    for iout in range(0,Nout) :
      isim=0
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        yt_funcs.plot_proj_gastemp(simdir, outnum, simtyp, boxlen_comov, overwrite)
        isim = isim + 1
  if(prop=='metalicity_proj') :
    boxlen_comov=6.77
    overwrite=False
    for iout in range(0,Nout) :
      isim=0
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        yt_funcs.plot_proj_gastemp(simdir, outnum, simtyp, boxlen_comov, overwrite)
        isim = isim + 1
  if(prop=='phaseplot') :
    boxlen_comov=6.77
    overwrite=False
    for iout in range(0,Nout) :
      isim=0
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        yt_funcs.tn_phaseplot(simdir, outnum, simtyp, boxlen_comov, overwrite)
        isim = isim + 1
  if(prop=='globalsink_vs_z')  :
    boxlen_comov = 6.77
    cfuncs.compare_sinkmasstot(simdir_arr, simtyp_arr, cbar_arr, linestyle_arr, boxlen_comov)
  if(prop=='haloprops') :
    overwrite_sink = False
    overwrite_star = False
    plotsinks      = True
    boxlen_comov = 6.77
   
    for iout in range(0,Nout) :
      isim = 0
      haloprop_outnum_arr = []
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout] 
        haloprop_outnum_arr.append(outnum)
        isim = isim + 1
      print('haloprop_outnum_arr :', haloprop_outnum_arr) 
      cfuncs.haloprops(simdir_arr, simtyp_arr, haloprop_outnum_arr, cbar_arr, linestyle_arr, 2000, boxlen_comov, overwrite_star, overwrite_sink, plotsinks)


  




   








