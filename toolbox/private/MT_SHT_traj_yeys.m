function ye_ys  = MT_SHT_traj_yeys(ye_ys_jj,l,N)
%PURPOSE  Reorder one trajectory sensitivity block for SHT boundary equations.
%
% The integrator stores d(ye)/d(ys) in trajectory ordering. This utility
% reshapes and reorders that block so it matches the segment and state ordering
% used by the MF- and HT-SHT Jacobian assembly.
%
% On input:
%   ye_ys_jj : Terminal-state sensitivity block for one segment.
%   l        : Number of sample points in the sensitivity block.
%   N        : First-order state dimension.
%
% On output:
%   ye_ys    : Reordered sensitivity block.
% Each trajectory couples only to its own initial state.
idx = (1:l)' + l*(0:N-1);
rows = repmat(idx,1,N);
cols = repelem(idx,1,N);
ye_ys = sparse(rows(:),cols(:),ye_ys_jj(:),l*N,l*N);
end
