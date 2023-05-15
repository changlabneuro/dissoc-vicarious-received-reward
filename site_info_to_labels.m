function s_info = site_info_to_labels(s_info)

dst_fields = {'lfp-channel', 'days', 'lfp-region', 'unit_index', 'spk-channel', 'spk-region'};

s_info(:, 3) = cellfun( @(x) sprintf('lfp-%s', x), s_info(:, 3), 'un', 0 );
s_info(:, 6) = cellfun( @(x) sprintf('spk-%s', x), s_info(:, 6), 'un', 0 );
s_info = fcat.from( s_info, dst_fields );

end