#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt

def main():
    pars_df = pd.read_csv(pars_path)
    site_df = pd.read_csv(site_path)

    big_tab = pd.merge( pars_df, site_df, left_on="Site", right_on="Label" )
    sub_tab = big_tab.query('Sampling=="in" and k==2')
    print sub_tab

    fig, (ax1,ax2,ax3) = plt.subplots(1, 3, sharey=True, figsize=(12,4))
    ax1.scatter(sub_tab["MAP"], sub_tab["Value"])
    ax2.scatter(sub_tab["MAT"], sub_tab["Value"])
    ax3.scatter(sub_tab["Latitude"], sub_tab["Value"])
    ax1.set_ylabel('k')
    ax1.set_xlabel('MAP (mm)')
    ax2.set_xlabel('MAT ($\degree$C)')
    ax3.set_xlabel('Latitude')
    plt.show()

    return None

if __name__=="__main__":
    pars_path = "../outputs/spring_parameters.csv"
    site_path = "../data/site_char.csv"

    main()
