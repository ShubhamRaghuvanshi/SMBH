import cfuncs
   
simdir_1='/u/sraghu/sim/SHUB/b256/sn_eq'
simdir_2='/u/sraghu/sim/SHUB/b256/sn_neq'
simdir_3='/u/sraghu/sim/SHUB/b256/sn_neq_18'
simdir_4='/u/sraghu/sim/SHUB/b256/sn_neq_lightstar'

simtyp_arr    = ['sn_neq_highres', 'sn_neq_lightstar']

simdir_arr  = [ simdir_3, simdir_4]
cbar_arr    = ['g', 'b']

#cfuncs.plot_mfblog(simdir_arr, simtyp_arr, cbar_arr, False)


cfuncs.plot_sfrlog(simdir_arr, simtyp_arr, cbar_arr, False)

#halonum = 1
#outnum_arr = [8,9]
#cfuncs.plot_halohist(simdir_arr, simtyp_arr, outnum_arr, cbar_arr, halonum)
