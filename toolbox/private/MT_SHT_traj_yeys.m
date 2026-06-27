function ye_ys  = MT_SHT_traj_yeys(ye_ys_jj,l,N)
ye_ys = zeros(l*N,l*N);
for kk = 1:l
    ye_ys(kk:l:end,kk:l:end) = reshape(ye_ys_jj(kk,:),N,N);
end
end