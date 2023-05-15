function [coh, labels, src_to_dst] = linearize_sfcoh(averageMat)

p = @(x) permute(x, [3, 1, 2]);

t1_cued_normSelf_popCoh     = p(squeeze(averageMat(:,:,1,1,:)));
t1_cued_normOther_popCoh    = p(squeeze(averageMat(:,:,1,2,:)));

t1_choice_normSelf_popCoh    = p(squeeze(averageMat(:,:,1,3,:)));
t1_choice_normOther_popCoh   = p(squeeze(averageMat(:,:,1,4,:)));

t2_cued_normSelf_popCoh     = p(squeeze(averageMat(:,:,2,1,:)));
t2_cued_normOther_popCoh    = p(squeeze(averageMat(:,:,2,2,:)));

t2_choice_normSelf_popCoh    = p(squeeze(averageMat(:,:,2,3,:)));
t2_choice_normOther_popCoh   = p(squeeze(averageMat(:,:,2,4,:)));

dir1 = { 'acc_bla' };
dir2 = { 'bla_acc' };
sl = { 'self-bottle' };
ol = { 'other-bottle' };
cued = { 'cued' };
choice = { 'choice' };

t1_cued_s = repmat( [dir1, sl, cued], size(t1_cued_normSelf_popCoh, 1), 1 );
t1_cued_o = repmat( [dir1, ol, cued], size(t1_cued_normOther_popCoh, 1), 1 );

t1_choice_s = repmat( [dir1, sl, choice], size(t1_choice_normSelf_popCoh, 1), 1 );
t1_choice_o = repmat( [dir1, ol, choice], size(t1_choice_normOther_popCoh, 1), 1 );

t2_cued_s = repmat( [dir2, sl, cued], size(t2_cued_normSelf_popCoh, 1), 1 );
t2_cued_o = repmat( [dir2, ol, cued], size(t2_cued_normOther_popCoh, 1), 1 );

t2_choice_s = repmat( [dir2, sl, choice], size(t2_choice_normSelf_popCoh, 1), 1 );
t2_choice_o = repmat( [dir2, ol, choice], size(t2_choice_normOther_popCoh, 1), 1 );

coh = [ 
  t1_cued_normSelf_popCoh;
  t1_cued_normOther_popCoh;
  
  t1_choice_normSelf_popCoh;
  t1_choice_normOther_popCoh;

  t2_cued_normSelf_popCoh;
  t2_cued_normOther_popCoh;

  t2_choice_normSelf_popCoh;
  t2_choice_normOther_popCoh;
];

labels = [
  t1_cued_s;
  t1_cued_o;

  t1_choice_s;
  t1_choice_o;

  t2_cued_s;
  t2_cued_o;

  t2_choice_s;
  t2_choice_o;
];
labels = fcat.from( labels, {'direction', 'outcome', 'trialtype'} );

src_to_dst = repmat( (1:size(averageMat, 5))', 8, 1 );
  
end