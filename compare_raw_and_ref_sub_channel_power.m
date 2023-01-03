cons = load( 'trial_data.mat' );
cons = cons.consolidated;

pl2_p = '/Volumes/external3/data/changlab/ptp-vicarious-reward/raw_pl2/Hitch_01262017.pl2';
pl2_f = shared_utils.io.filenames( pl2_p, true );

chan_info = cons.pl2_info.channel_map(pl2_f);
chans = arrayfun( @(x) x.channels(1), chan_info );
regs = arrayfun( @(x) x.region, chan_info, 'un', 0 );

m = [];
fs = [];
for i = 1:numel(regs)
  chan_str = dsp3.channel_n_to_str( 'FP', chans(i) );
  ad = PL2Ad( pl2_p, chan_str );
  if ( isempty(m) )
    m = nan( numel(regs), numel(ad.Values) );
  end
  m(i, :) = ad.Values;
  fs = ad.ADFreq;
end

%%

file_ind = strcmp( cons.pl2_info.files, pl2_f );
start_t = cons.pl2_info.start_times(file_ind);
t_series = (0:size(m, 2)-1) / fs + start_t;

event_labs = fcat.from( cons.events.labels );
align_labs = fcat.from( cons.align.labels );
event_ind = find( event_labs, cons.pl2_info.sessions(file_ind) );
align_ind = find( align_labs, cons.pl2_info.sessions(file_ind) );

fix_ts = cons.events.data(event_ind, cons.event_key('fixOn'));
fix_ts(fix_ts == 0) = nan;
align_plex = cons.align.data(align_ind, cons.align_key('plex'));
align_picto = cons.align.data(align_ind, cons.align_key('picto'));

plex_fix_ts = shared_utils.sync.cinterp( fix_ts, align_picto, align_plex );

%%

ref_ind = ismember( regs, 'ref' );
assert( nnz(ref_ind) == 1 );
rest_inds = setdiff( 1:numel(regs), find(ref_ind) );

ref_sub_m = m(rest_inds, :) - m(ref_ind, :);

bin_w = 150;
bin_s = 50;
nfft = 256;
f1 = 2.5;
f2 = 250;
chronux_params = shared_utils.general.get( dsp3.make.defaults.coherence, 'chronux_params' );
do_filt = true;

% power_func = @(x) periodogram( x, [], nfft, fs );
power_func = @(x) periodogram( x, [], 1:0.5:20, fs );
% power_func = @(x) mtspectrumc( x, chronux_params );

inds = shared_utils.vector.slidebin( 1:size(m, 2), bin_w, bin_s, true );
inds = inds((1:1e4) + 0);

psd_base = [];
psd_ref_sub = [];

binned_t_series = cellfun( @(x) min(t_series(x)), inds );

for i = 1:numel(inds)
  fprintf( '\n %0.3f', i/numel(inds) * 1e2 );
  
  ind = inds{i};
  
  for j = 1:numel(rest_inds)
    ri = rest_inds(j);
    
    m_win = m(ri, ind);
    ref_m_win = ref_sub_m(j, ind);
    
    if ( do_filt )
      m_win = dsp3.zpfilter( m_win, f1, f2, fs );
      ref_m_win = dsp3.zpfilter( ref_m_win, f1, f2, fs );
    end
    
    [one_psd_base, f] = power_func( m_win );
    one_psd_ref = power_func( ref_m_win );
    
    keep_f = find( f < 12 );
    if ( isempty(psd_base) )
      psd_base = nan( numel(rest_inds), numel(keep_f), numel(inds) );
      psd_ref_sub = nan( size(psd_base) );
    end
    
    psd_base(j, :, i) = one_psd_base(keep_f);
    psd_ref_sub(j, :, i) = one_psd_ref(keep_f);
    f = f(keep_f);
  end
end

ts = cellfun(@median, inds) / fs;

%%

nearest_bins = nan( numel(plex_fix_ts), 1 );
for i = 1:numel(plex_fix_ts)
  if ( isnan(plex_fix_ts(i)) || plex_fix_ts(i) > max(binned_t_series) ), continue; end
  [~, nearest_bins(i)] = min( abs(plex_fix_ts(i) - binned_t_series) );
end

nearest_bins = nearest_bins(~isnan(nearest_bins));
base_psd_base = nanmean( psd_base(:, :, nearest_bins), 3 );
base_psd_ref = nanmean( psd_ref_sub(:, :, nearest_bins), 3 );

psd_base_norm = psd_base ./ base_psd_base;
psd_ref_norm = psd_ref_sub ./ base_psd_ref;

%%

do_norm = true;

if ( do_norm )
  psd_plt = [ psd_base_norm; psd_ref_norm ];
else
  psd_plt = db( [psd_base; psd_ref_sub] );
end

labels = [ repmat(columnize(regs(rest_inds)), 2, 1) ...
  , {'raw', 'raw', 'ref-sub', 'ref-sub'}' ];

t_ind = binned_t_series > min(binned_t_series) + 300 & ...
  binned_t_series < max(binned_t_series) - 150;

figure(1);
axs = plots.panels( size(psd_plt, 1), true );
for i = 1:numel(axs)
  ax = axs(i);
  to_plt = squeeze( psd_plt(i, :, t_ind) );
  to_plt = imgaussfilt( to_plt, 5 );
  imagesc( ax, binned_t_series(t_ind), f, to_plt );
  set( ax, 'ydir', 'normal' );
  title( ax, strjoin(labels(i, :), ' | ') );
  colorbar( ax );
  ylabel( ax, 'Hz' );
  xlabel( ax, 'Time from session start (s)' );
end

%%

if ( do_norm )
  shared_utils.plot.set_clims( axs, [0, 20] );
else
  % shared_utils.plot.set_clims( axs, [0, 10e-3] );
  shared_utils.plot.set_clims( axs, [-200, -10] );
end
