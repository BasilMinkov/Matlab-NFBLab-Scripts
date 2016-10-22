function chanlocs_vis = makeChanlocsVis()
load('chanlocs_mod.mat')
used_ch = { 'Fp1','Fp2','F7','F3','Fz','F4','F8','Ft9','Fc5','Fc1','Fc2','Fc6','Ft10','T7','C3','Cz',...
'C4','Tp9','Cp5','Cp1','Cp2','Cp6','Tp10','P7','P3','Pz','P4','P8','O1','Oz','O2'};
[chan_idx_chanlocs,ifvis] = find_str_idx(lower({chanlocs.labels}), lower(used_ch));
chanlocs_vis = chanlocs(chan_idx_chanlocs);
end
