%% To extract the voltage from 32 channels and timestamps for stimuli for experiment 20211116

Record1_data_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys1_2021-11-11_11-41-13\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record1_stim_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys1_2021-11-11_11-41-13\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record2_data_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys2_2021-11-11_12-00-58\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record2_stim_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys2_2021-11-11_12-00-58\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record3_data_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys3_2021-11-11_12-31-48\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record3_stim_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys3_2021-11-11_12-31-48\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record4_data_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys4_2021-11-11_12-52-56\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record4_stim_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys4_2021-11-11_12-52-56\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record5_data_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys5_2021-11-11_13-12-55\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record5_stim_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys5_2021-11-11_13-12-55\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record6_data_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys6_2021-11-11_13-32-23\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record6_stim_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys6_2021-11-11_13-32-23\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record7_data_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys7_2021-11-11_13-51-46\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record7_stim_20211116 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211116_AAV8_hsyn_GFP\NT_Ephys7_2021-11-11_13-51-46\Record Node 110\experiment1\recording1\structure.oebin','events',1);

%% To save the voltage from 32 channels and timestamps for stimuli for experiment 20211021
save('Record1_20211116','Record1_data_20211116','Record1_stim_20211116','-v7.3');
save('Record2_20211116','Record2_data_20211116','Record2_stim_20211116','-v7.3');
save('Record3_20211116','Record3_data_20211116','Record3_stim_20211116','-v7.3');
save('Record4_20211116','Record4_data_20211116','Record4_stim_20211116','-v7.3');
save('Record5_20211116','Record5_data_20211116','Record5_stim_20211116','-v7.3');
save('Record6_20211116','Record6_data_20211116','Record6_stim_20211116','-v7.3');
save('Record7_20211116','Record7_data_20211116','Record7_stim_20211116','-v7.3');
