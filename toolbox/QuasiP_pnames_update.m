function data = QuasiP_pnames_update(tbid, data, pnames)
if strcmp(tbid, 'QuasiP')
    prefix = '';
else
    parts = strsplit(tbid, '.');
    prefix = [parts{1}, '.'];
end
data.pnames_np  = pnames;
data.my_tbid    = tbid;
data.prefix     = prefix;
%----------------------------------------------------------
if ~isempty(pnames) && iscell(pnames)
    for i = 1:length(pnames)
        data.pnames{i} = [prefix, pnames{i}];
    end
else
    data.pnames = pnames;
end
%----------------------------------------------------------
data.kOMnames = cell(1, data.ode_opts.e);
for ii = 1 : data.ode_opts.e
    data.kOMnames{ii} = sprintf([prefix 'k_OM%d'], ii);
end
%----------------------------------------------------------
nUOM = data.ode_opts.d - data.ode_opts.e;
data.uOMnames = cell(1, nUOM);
for ii = 1 : nUOM
    data.uOMnames{ii} = sprintf([prefix 'u_om%d'], data.ode_opts.e + ii);
end
end
