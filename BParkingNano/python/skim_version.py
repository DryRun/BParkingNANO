# skim_version = "v1_0" 
# skim_version = "v1_1" # Redo Bs skim with fixed MC sample
# skim_version = "v1_2" # Fix gen particles, 531 was not in the list
#skim_version = "v1_3" # Add all leptons and their parents+grandparents to the saved gen particles
#skim_version = "v1_4" # Lower pT cut on Bs Ks to 0.5, symmetrize phi window to [0.96, 1.08]. Also, Bs to KKuu probefilter has more stats.
#skim_version = "v1_5" # Loosed pre-vtx phi window to [0.52, 1.52]. Set post-vtx phi window to [0.96, 1.08].
#skim_version = "v1_5_veryloose" # Remove almost all vtx reqs to debug efficiency loss
#skim_version = "v1_6" # More pre-vtx loosening, e.g. trk eta < 2.5, not 2.4
#skim_version = "v1_6_lzma6" # Test compression level 6
#skim_version = "v1_6_zlib6" # Test ZLIB compression level 6
#skim_version = "v1_7" # v1_6 data had 200 LS/job, which caused timeouts. Try 50 instead (might have to hadd)
#skim_version = "v1_7_1" # Redo 1_7 failed jobs, which expired
#skim_version = "v1_7_2" # Bug-fix and ext Bs_probefilter MC

#skim_version="v2_0" # Add Jpsi mass constraint.
#skim_version="v2_1" # Fix ordering of tracks passed to KinVtxFitter.
skim_version="v2_2" # Fix output data - leptons and tracks were mixed up.