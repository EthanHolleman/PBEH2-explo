from pathlib import Path

LOW_DIR = '/group/flchedingrp/mitochi/Work/Project/Ethan/240325/5_Analysis/2_deleted/67_72_PFC8TACINIT_240406/PBEH2_240409_L1_t45_w15_l50/PEAKS_GENOME'
HIGH_DIR = '/group/flchedingrp/mitochi/Work/Project/Ethan/240325/5_Analysis/2_deleted/67_72_PFC8TACINIT_240406/PBEH2_240409_L1_t55_w20_l100/PEAKS_GENOME'

low_files = list(Path(LOW_DIR).iterdir())
high_files = list(Path(HIGH_DIR).iterdir())

def just_names(file_list):

    return [f.name for f in file_list]


low_names = set(just_names(low_files))
high_names = set(just_names(high_files))

shared = low_names.intersection(high_names)
diff = low_names.difference(high_names)

print(len(shared))