import os


# if laptop: "/home/franzese/projects/hotspot_signature_panel/data/"
# if workstation: "/fs/cbcb-lab/mdml/users/franzese/projects/signature-panel/signature-panel/data/"

GLOBAL_DATA_DIR = "/home/franzese/projects/hotspot_signature_panel/data/"
GLOBAL_CHR_MTX_DIR = GLOBAL_DATA_DIR + "individual_chromosome_matrices/"

GLOBAL_BASELINE_SBS_DIR = GLOBAL_DATA_DIR + "BASELINE_PANELS/SBS_MATRICES/"
GLOBAL_BASELINE_EST_DIR = GLOBAL_DATA_DIR + "BASELINE_PANELS/SIGNATURE_ESTIMATES/"


GLOBAL_OUT_DIR = "/home/franzese/projects/hotspot_signature_panel/out/" 

GLOBAL_PANEL_SBS_DIR = GLOBAL_OUT_DIR + "SIGNATURE_PANELS/SBS_MATRICES/"
GLOBAL_PANEL_EST_DIR = GLOBAL_OUT_DIR + "SIGNATURE_PANELS/SIGNATURE_ESTIMATES/"


def estimate_signatures_for_panels(baseline_flag=True):
    
    if baseline_flag:
        sbs_dir = GLOBAL_BASELINE_SBS_DIR
        est_dir = GLOBAL_BASELINE_EST_DIR
    else:
        sbs_dir = GLOBAL_PANEL_SBS_DIR
        est_dir = GLOBAL_PANEL_EST_DIR

    sbs_file_ls = os.listdir(sbs_dir)
    

    for f in sbs_file_ls:
        
        outfile = est_dir + f[:-4] + "_SIG_EST.tsv"
        os.system("python signature-estimation-py/signature_estimation.py -mf " + 
                sbs_dir + f + " -sf data/cosmic-signatures.tsv -of " + outfile)

    return

estimate_signatures_for_panels(False)

