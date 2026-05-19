%% To extract the voltage from 32 channels and timestamps for stimuli for experiment 20211110

Record1_data_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys1_2021-11-10_12-06-16\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record1_stim_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys1_2021-11-10_12-06-16\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record2_data_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys2_2021-11-10_12-25-51\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record2_stim_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys2_2021-11-10_12-25-51\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record3_data_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys3_2021-11-10_12-45-12\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record3_stim_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys3_2021-11-10_12-45-12\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record4_data_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys4_2021-11-10_13-04-16\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record4_stim_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys4_2021-11-10_13-04-16\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record5_data_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys5_2021-11-10_13-23-35\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record5_stim_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys5_2021-11-10_13-23-35\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record6_data_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys6_2021-11-10_13-43-07\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record6_stim_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys6_2021-11-10_13-43-07\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record7_data_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys7_2021-11-10_14-02-20\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record7_stim_20211110 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\20211110_AAV8_hsyn_ChR2_x6\NT_ephys7_2021-11-10_14-02-20\Record Node 110\experiment1\recording1\structure.oebin','events',1);

%% To save the voltage from 32 channels and timestamps for stimuli for experiment 20211021
save('Record1_20211110','Record1_data_20211110','Record1_stim_20211110','-v7.3');
save('Record2_20211110','Record2_data_20211110','Record2_stim_20211110','-v7.3');
save('Record3_20211110','Record3_data_20211110','Record3_stim_20211110','-v7.3');
save('Record4_20211110','Record4_data_20211110','Record4_stim_20211110','-v7.3');
save('Record5_20211110','Record5_data_20211110','Record5_stim_20211110','-v7.3');
save('Record6_20211110','Record6_data_20211110','Record6_stim_20211110','-v7.3');
save('Record7_20211110','Record7_data_20211110','Record7_stim_20211110','-v7.3');
