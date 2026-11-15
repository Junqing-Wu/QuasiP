function QuasiP_print_equation_table(data)
%QUASIP_PRINT_EQUATION_TABLE  Print the assembled model and PHI summary.
%
% Build a readable string for the selected first- or second-order equation and
% print the frequency counts, parameter names, period factors, mesh sizes, and
% discretization dimensions. This provides a reproducible record of the
% problem that was passed to the private residual.
%
% On input:
%   data : Initialized QuasiP data structure.
%
% On output:
%   None. The equation table is printed to the MATLAB command window.
if data.isL
    terms = [];names = {'A', 'B', 'C'};exps = {'*ddx', '*dx', '*x'};
    for i = 1:3
        if data.isABC(i)
            term = names{i};
            if data.isABCwithp(i)
                term = [term '(p)'];
            end
            term = [term exps{i}];
            if isempty(terms)
                terms = term;
            else
                terms = [terms ' + ' term];
            end
        end
    end
    terms = [terms ' = f(x'];
else
    terms = 'dx = f(x';
end
%----------------------------------------------------------
if data.isfwithu
    terms = [terms ',dx'];
end
if data.isfwithp
    terms = [terms ',p'];
end
for ii = 1:data.ode_opts.e
    terms = [terms, sprintf(',th%d', ii)];
end
terms = [terms, ')'];
%----------------------------------------------------------
headers = {'d','e','pnames','kOMs','uoms','PHI','P','N','L','P*S'};
switch data.PHI.param
    case{'MP'}
        values  = {
            num2str(data.ode_opts.d), ...
            num2str(data.ode_opts.e), ...
            strjoin(data.pnames, '+'), ...
            strjoin(data.kOMnames, '+'), ...
            strjoin(data.uOMnames, '+'), ...
            strcat(string(data.PHI.param), '_',string(data.PHI.method), '_', strjoin(data.PHI.type,'+'), '_Od(', strjoin(string(data.PHI.OmOd(:)),','), ')'), ...
            strjoin(string(data.PHI.P), '+'), ...
            num2str(data.ode_opts.n), ...
            num2str(data.PHI.L), ...
            strjoin(string(data.PHI.pS), '+')
            };
    case{'MT'}
        switch data.PHI.method
            case{'DM'}

            case{'SHT'}
                pS = cell(1,data.ode_opts.d);
                if data.PHI.IsFixedStep
                    pS{data.PHI.Omt}   = data.PHI.P{data.PHI.Omt}*data.PHI.L1;
                    for ii = 1:data.ode_opts.d-1
                        pS{data.PHI.Omnot(ii)} = data.PHI.pSnot{ii};
                    end
                else
                    pS{data.PHI.Omt}   = 'UN';
                    for ii = 1:data.ode_opts.d-1
                        pS{data.PHI.Omnot(ii)} = data.PHI.pSnot{ii};
                    end
                end
                values  = {
                    num2str(data.ode_opts.d), ...
                    num2str(data.ode_opts.e), ...
                    strjoin(data.pnames, '+'), ...
                    strjoin(data.kOMnames, '+'), ...
                    strjoin(data.uOMnames, '+'), ...
                    [data.PHI.param '_' data.PHI.method '_' data.PHI.type '_' func2str(data.PHI.ODEsolver) '_t(' num2str(data.PHI.Omt) ')_M(' num2str(data.PHI.M) ')'], ...
                    strjoin(string(data.PHI.P), '+'), ...
                    num2str(data.PHI.N), ...
                    num2str(data.PHI.L), ...
                    strjoin(string(pS), '+')
                    };
        end
end
%----------------------------------------------------------
values_char = cellfun(@(v) char(v), values, 'UniformOutput', false);
col_widths = cellfun(@(h,v) max(strlength(h), strlength(v))+4, headers, values_char);
col_widths(end) = col_widths(end)-4;
total_width = sum(col_widths);
line = repmat('=', 1, total_width);
term_str = [' ' terms ' '];
fill_len = total_width  - length(term_str);
left_fill = repmat('>', 1, ceil(fill_len/2));
right_fill = repmat('<', 1, floor(fill_len/2));
title_line = [left_fill term_str right_fill];
%----------------------------------------------------------
disp(line);
disp(title_line);
disp(line);
fmt = strjoin(compose('%%-%ds', col_widths), '');
fprintf([fmt '\n'], headers{:});
fprintf([fmt '\n'], values{:});
disp(line);
end