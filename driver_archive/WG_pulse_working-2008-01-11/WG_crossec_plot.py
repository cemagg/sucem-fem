line_types = ('k-', 'k:', 'k-.', 'k--')

res_2b = pickle.load(file('+coupledb_wg_direch_2nd_crossection.pickle'))[2][32]
res_3b = pickle.load(file('+coupledb_wg_direch_3rd_crossection.pickle'))[3][32]
res_2a = pickle.load(file('+coupleda_wg_direch_2nd_crossection.pickle'))[2][32]
an_res = pickle.load(file('+WG_analytical_result_9s_div32_new.pickle'))
dt = analytical_dt = an_res.dt
cf_2b = res_2b.crosssec_field
cf_3b = res_3b.crosssec_field
cf_2a = res_2a.crosssec_field

peak=dict()
peak['2a'] = 973
peak['2b'] = 983
peak['an'] = 983


dx = 1/100.
x = N.arange(101)*dx
cf_an = N.sin(N.pi*x)
# plot(x, cf_an*an_res.ts[peak['an']])
# plot(x, cf_2b[peak['2b'],:,1])
# plot(x, cf_2a[peak['2a'],:,1])

t = N.arange(len(an_res.ts))*dt
r_s, r_e = int(1/dt), int(3.5/dt)
lt = line_types.__iter__()
plot(t[r_s:r_e], res_2b.ts_modeintg_n[r_s:r_e], lt.next(),label='EB Mode Integrated')
plot(t[r_s:r_e], res_2b.point_reconstructed_ts[r_s:r_e], lt.next(), label='EB Field Point')
plot(t[r_s:r_e], res_2a.ts_modeintg_n[r_s:r_e], lt.next(), label='EBHD Mode Integrated')
plot(t[r_s:r_e], res_2a.point_reconstructed_ts[r_s:r_e], lt.next(), label='EBHD Field Point')
xlabel('Time (s)')
ylabel('E field y-component (V/m)')
legend()
