is_cluster = true;
dir_method = 'pdc';
do_site_mean = true;
do_n_match = true;
allow_overwrite = true;

max_model_order = 40;
% max_model_order = 10;
% max_model_order = 40;

% max_num_files = inf;
max_num_files = 5;

verbose = ~is_cluster;

if ( is_cluster )  
  rd = '/home/naf3/Repositories';
  run( fullfile(rd, 'asympPDC', 'startup.m') );

  parpool( feature('numcores') );

  addpath( fullfile(rd, 'ptp-vicarious-reward') );
  addpath( genpath(fullfile(rd, 'dsp')) );

  dr = '/gpfs/milgram/project/chang/CHANG_LAB/ptp8/rewardLFP';
  dst_p = '/gpfs/milgram/project/chang/CHANG_LAB/naf3/Data/ptp-vicarious-reward';

  lfp_mats = shared_utils.io.findmat( fullfile(dr, 'reward-aligned') );
else
  repadd( 'ptp-vicarious-reward', true );

  dst_p = '/Volumes/external3/data/changlab/ptp-vicarious-reward';
  lfp_mats = shared_utils.io.findmat( fullfile(dst_p, 'reward_aligned_lfp', 'reward-aligned') );
end

switch ( dir_method )
  case 'pdc'
    pdc_subdir = 'pdc';
  case 'pdc-orig'
    pdc_subdir = 'pdc-orig';
  case 'dtf'
    pdc_subdir = 'dtf';
  otherwise
    error( 'Unrecognized method "%s".', dir_method );
end

sm_str = '';
if ( do_site_mean ), sm_str = '-site-mean'; end

nm_str = '';
if ( do_n_match ), nm_str = '-nmatch'; end

is_orig_pdc = strcmp( dir_method, 'pdc-orig' );

pdc_subdir = sprintf( '%s%s%s-mo-%d', pdc_subdir, sm_str, nm_str, max_model_order );
pdc_dst_p = fullfile( dst_p, pdc_subdir );
shared_utils.io.require_dir( pdc_dst_p );

site_pairs = shared_utils.io.fload( fullfile(dst_p, 'site-pairs/pairs.mat') );
site_mean_each = {'sites', 'channels', 'regions', 'outcomes', 'trialtypes', 'administration'};
match_n_within = { 'sites', 'channels', 'regions' };
match_n_across = { 'outcomes' };

%%

parfor m = 1:min(max_num_files, numel(lfp_mats))
  
if ( is_cluster )
  warning( 'off', 'all' );
end

%%

fprintf( '\n %d of %d', m, numel(lfp_mats) );

%%

lfp_file = load( lfp_mats{m} );
  
%%

[data, labels, categories, t] = signal_container_lfp_to_unformatted( lfp_file.lfp );
labels = fcat.from( labels, categories );

mask = pipe( rowmask(labels) ...
  , @(m) findnone(labels, 'errors', m) ...
  , @(m) find(labels, 'pre', m) ...
);

if ( do_n_match )
  match_I = findall( labels, setdiff(site_mean_each, union(match_n_across, match_n_within)), mask );
  mask = n_match( labels, match_I, match_n_within, match_n_across );
end

[pair_sets, pair_labels] = find_pairs( labels, 'regions', 'channels', mask );

day_str = combs( labels, 'days' );
assert( numel(day_str) == 1 );
target_pairs = site_pairs.channels{strcmp(site_pairs.days, day_str)};

%%

bi = shared_utils.vector.slidebin( 1:size(data, 2), 150, 50, true );

for pi = 1:size(pair_sets, 1)
  fprintf( '\n\t %d of %d', pi, size(pair_sets, 1) );
  target_pair = pair_sets(pi, :);
  
  if ( ~accept_pair(target_pairs, pair_labels(pi, :)) )
    continue; 
  end
  
  dst_fp = fullfile( pdc_dst_p, sprintf('pdc-%d-%d.mat', m, pi) );
  if ( ~allow_overwrite && exist(dst_fp, 'file') )
    continue;
  end
  
  if ( 0 )  % debug: use only a couple trials
    target_pair = eachcell( @(x) x(1:2), target_pair );
  end

  [pdc, pdc_ls, pdc_f] = run_pdc_pair( data, labels, target_pair, bi ...
    , 'verbose', verbose ...
    , 'asymp_pdc_params', {'dir_method', dir_method, 'maxIP', max_model_order} ...
    , 'is_orig_pdc_method', is_orig_pdc ...
  );
  pdc_t = cellfun( @(x) t(x(1)), bi );
  
  if ( do_site_mean )
    [pdc_ls, I] = keepeach( pdc_ls, site_mean_each );
    pdc = bfw.row_nanmean( pdc, I );
  end
  
  do_save( dst_fp, pdc, pdc_ls, pdc_f, pdc_t );
end

end

%%

function tf = accept_pair(target_pairs, current_pair)

sp = sortrows( target_pairs' )';
cp = sortrows( current_pair(:) )';
tf = ismember( string(cp), string(sp), 'rows' );

end

function do_save(p, pdc, pdc_ls, pdc_f, pdc_t)
save( p, 'pdc', 'pdc_ls', 'pdc_f', 'pdc_t' );
end