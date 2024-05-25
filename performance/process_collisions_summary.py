import pandas as pd 
import numpy as np
from scipy.stats import gmean
import itertools
import math


df = pd.read_csv("collisions_summary.csv")
# Remove useless experiments! 
df = df[df['C_size'] > df['log2(nwords)']]

# Get all the possible values of C
C_sizes = df["C_size"].unique()
C_sizes = sorted(C_sizes) # sort them for better visuals

# print(f"C_sizes are {C_sizes}")
# print(f"Column Names are {df.columns.values.tolist()}")

# Ignore experiments where we have more memory than the problem size
df = df[df['C_size'] > df['log2(nwords)']]
# because in our experimentatino we use a multiple power of two 
df['log2(nwords)'] = df['log2(nwords)'].astype(int)
print(f"Raw data has {df.shape[0]} samples")

# Get all C_sizes 
C_sizes = df["C_size"].unique()
nwords = df["log2(nwords)"].unique()
nwords = sorted(nwords)


# npoints
# factor of sqrt(n^3/w)
# sqrt(n^3/w) (theory)
# standard deviation 
# sample size or nruns
# nfunction updates 
# npoints / ndist 


# main row 
row_nwords = {"nwords": ["2^" + str(w) for w in nwords]}

factor_n3w = {}
factor_n3w.update(row_nwords)

stds = {}
stds.update(row_nwords)

avgs = {}
avgs.update(row_nwords)


nsamples = {}
nsamples.update(row_nwords)

nfunction_updates = {}
nfunction_updates.update(row_nwords)

# npts_over_ndist = {}
# npts_over_ndist.update(row_nwords)


def dict_to_row(d, name, indices):
    row = {name: []}
    for i in indices:
        try: # if value exist
            row[name].append(d[i])
        except: # add -1 to indicate no value
            row[name].append(-1)
    return row

# print("avg stats for a normal run...\n"
#      "scale: sqrt(n^3/mem)\n"
#      "----------------------\n")


for c in C_sizes:
    # print(f"*** C_size = {c}, scale = 2^{c*3/2:.2f}/sqrt(mem) ***\n")
    df_C = df[df["C_size"] == c]
    print(f"N=2^{c} has {df_C.shape[0]} samples")
    difficulties = df_C["difficulty"].unique()
    difficulties = sorted(difficulties)
    
    rams = df_C["log2(nwords)"].unique()
    rams = sorted(rams)

    row_factor_n3w = dict({})
    row_stds = dict({})
    row_avgs = dict({})
    row_nsamples = dict({})
    row_nfunction_updates = dict({})     
    # row_npts_over_ndist = dict({})
    
    for d in difficulties[::-1]: # better read them in reverse order 
        df_C_d = df_C[df_C['difficulty'] == d]
        #print(f"******* has {df_C_d.shape[0]} samples. diff={d}")
        for m in rams:
            df_C_d_m = df_C_d[df_C_d['log2(nwords)'] == m]
            if df_C_d_m.empty: # skip this dataframe since we are combining 
                continue # all values using brute force
            print(f"---------> {df_C_d_m.shape[0]} samples, nwords=2^{m}, diff={d}")
 
                
            # Calculate the mean of the desired column (e.g., 'Column1')
            # ndist_avg = gmean(df_C_d_m["#distinguished_points"])
            npoints_avg = gmean(df_C_d_m["#points"])
            # collisions_avg = gmean(df_C_d_m["#collisions"])
            
            # in the papaer they have divided by sqrt(n^3/w)
            lg2_scale_factor = (c*3 - m)/2
            scale_factor = 2**lg2_scale_factor
            
            # ndist_scaled = ndist_avg / scale_factor
            npoints_scaled =  npoints_avg / scale_factor
            # ncollisions_scaled = collisions_avg / scale_factor
        
            # ndist_std_scaled = df_C_d_m["#distinguished_points"].std() / scale_factor
            log2_npoints = np.log2(df_C_d_m["#points"]) # row
            log2_npoints_std =  log2_npoints.std() # single value
            
            # npoints_std_scaled = (2**log2_npoints_std) / scale_factor
            # ncollisions_std_scaled = df_C_d_m["#collisions"].std() / scale_factor

            nupdates_avg = df_C_d_m["#updates"].mean()
            # add it to a table to be saved in tables.txt
            row_factor_n3w["2^" + str(m) ] = npoints_scaled
            row_stds["2^" + str(m) ]       = log2_npoints_std
            row_avgs["2^" + str(m) ]       = f"2^{np.log2(npoints_avg):.2f}"
            row_nsamples["2^" + str(m) ]   = df_C_d_m.shape[0]
            row_nfunction_updates["2^" + str(m) ] = f"2^{np.log2(nupdates_avg):.2f}"
            # row_npts_over_ndist["2^" + str(m) ]   = npoints_avg / ndist_avg
    
    factor_n3w.update(dict_to_row(row_factor_n3w, f"2^{c}", row_nwords["nwords"]))
    stds.update(dict_to_row(row_stds, f"2^{c}", row_nwords["nwords"]))
    avgs.update(dict_to_row(row_avgs, f"2^{c}", row_nwords["nwords"]))
    nsamples.update(dict_to_row(row_nsamples, f"2^{c}", row_nwords["nwords"]))

    nfunction_updates.update(dict_to_row(row_nfunction_updates, f"2^{c}", row_nwords["nwords"]))
    # npts_over_ndist.update(dict_to_row(row_npts_over_ndist, f"2^{c}", row_nwords["nwords"]))
    
    # print("="*40 + "\n")