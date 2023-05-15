function [pdcResults] = F2F_PDC_Lite_func(input_labels,input_data)

offsetMs = 1000; 
pdc_timeRange = [-600,600];

choice_idxs = find(input_labels, 'choice');
cued_idxs = find(input_labels, 'cued');
acc_idxs = find(input_labels, 'acc');
bla_idxs = find(input_labels, 'bla');

self_idxs = find(input_labels, 'self');
other_idxs = find(input_labels, 'other');
both_idxs = find(input_labels, 'both');
none_idxs = find(input_labels, 'none');

choice_both_idxs = intersect(both_idxs, choice_idxs);
choice_self_idxs = intersect(self_idxs, choice_idxs);
choice_other_idxs = intersect(other_idxs, choice_idxs);
choice_none_idxs = intersect(none_idxs, choice_idxs);

cued_both_idxs = intersect(both_idxs, cued_idxs);
cued_self_idxs = intersect(self_idxs, cued_idxs);
cued_other_idxs = intersect(other_idxs, cued_idxs);
cued_none_idxs = intersect(none_idxs, cued_idxs);

bla_choice_both_idxs = intersect(choice_both_idxs, bla_idxs);
bla_choice_self_idxs = intersect(choice_self_idxs, bla_idxs);
bla_choice_other_idxs = intersect(choice_other_idxs, bla_idxs);
bla_choice_none_idxs = intersect(choice_none_idxs, bla_idxs);

bla_cued_both_idxs = intersect(cued_both_idxs, bla_idxs);
bla_cued_self_idxs = intersect(cued_self_idxs, bla_idxs);
bla_cued_other_idxs = intersect(cued_other_idxs, bla_idxs);
bla_cued_none_idxs = intersect(cued_none_idxs, bla_idxs);

acc_choice_both_idxs = intersect(choice_both_idxs, acc_idxs);
acc_choice_self_idxs = intersect(choice_self_idxs, acc_idxs);
acc_choice_other_idxs = intersect(choice_other_idxs, acc_idxs);
acc_choice_none_idxs = intersect(choice_none_idxs, acc_idxs);

acc_cued_both_idxs = intersect(cued_both_idxs, acc_idxs);
acc_cued_self_idxs = intersect(cued_self_idxs, acc_idxs);
acc_cued_other_idxs = intersect(cued_other_idxs, acc_idxs);
acc_cued_none_idxs = intersect(cued_none_idxs, acc_idxs);

pdcResults.sessionID = input_labels('session_ids');
pdcResults.monkey = input_labels('monkeys');
pdcResults.monkey = input_labels('days');
pdcResults.chan1Region  = 'BLA';
pdcResults.chan2Region  = 'ACC';
pdcResults.timeRangeMs = pdc_timeRange;


%% Choice trials
choice_both_data(:,:, 1) = input_data(bla_choice_both_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
choice_both_data(:,:, 2) = input_data(acc_choice_both_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
[~, pdcResults.choice_both_PDCavg] = partialDirectedCOH(choice_both_data);

choice_self_data(:,:, 1) = input_data(bla_choice_self_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
choice_self_data(:,:, 2) = input_data(acc_choice_self_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
[~, pdcResults.choice_self_PDCavg] = partialDirectedCOH(choice_self_data);

choice_other_data(:,:, 1) = input_data(bla_choice_other_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
choice_other_data(:,:, 2) = input_data(acc_choice_other_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
[~, pdcResults.choice_other_PDCavg] = partialDirectedCOH(choice_other_data);

choice_none_data(:,:, 1) = input_data(bla_choice_none_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
choice_none_data(:,:, 2) = input_data(acc_choice_none_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
[~, pdcResults.choice_none_PDCavg] = partialDirectedCOH(choice_none_data);


%% Cued
cued_both_data(:,:, 1) = input_data(bla_cued_both_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
cued_both_data(:,:, 2) = input_data(acc_cued_both_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
[~, pdcResults.cued_both_PDCavg] = partialDirectedCOH(cued_both_data);

cued_self_data(:,:, 1) = input_data(bla_cued_self_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
cued_self_data(:,:, 2) = input_data(acc_cued_self_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
[~, pdcResults.cued_self_PDCavg] = partialDirectedCOH(cued_self_data);

cued_other_data(:,:, 1) = input_data(bla_cued_other_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
cued_other_data(:,:, 2) = input_data(acc_cued_other_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
[~, pdcResults.cued_other_PDCavg] = partialDirectedCOH(cued_other_data);

cued_none_data(:,:, 1) = input_data(bla_cued_none_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
cued_none_data(:,:, 2) = input_data(acc_cued_none_idxs, pdc_timeRange(1)+offsetMs : pdc_timeRange(2)+offsetMs);
[~, pdcResults.cued_none_PDCavg] = partialDirectedCOH(cued_none_data);

