function time_freq_axis_labels(axs, t, f, t_spc, f_spc)

if ( nargin < 5 || isepty(f_spc) )
  f_spc = 2;
end
if ( nargin < 4 || isempty(t_spc) )
  t_spc = 4;
end

shared_utils.plot.fseries_yticks( axs, round(flip(f)), f_spc );
shared_utils.plot.tseries_xticks( axs, round(t), t_spc );

end