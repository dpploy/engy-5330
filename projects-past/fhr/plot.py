#!/usr/bin/env python
import matplotlib.pyplot as plt
import pandas as pd

def main():
    #df1 = pd.read_csv('reaction_1D_steady_out_omega_1_0002.csv')
    #df2 = pd.read_csv('reaction_1D_steady_out_omega_2_0002.csv')
    #df1 = pd.read_csv('reaction_1D_steady-added-source_out_omega_1_0002.csv')
    #df2 = pd.read_csv('reaction_1D_steady-added-source_out_omega_2_0002.csv')
    df1 = pd.read_csv('heat-modified_out_omega_1_0002.csv')
    df2 = pd.read_csv('heat-modified_out_omega_2_0002.csv')

    (fig, ax1) = plt.subplots(1, figsize=(8, 4))

    ax1.plot(df1['x'], df1['u'],'r*-')
    ax1.plot(df2['x'], df2['v'],'r*--')

    plt.grid(True)
    plt.savefig('plot.png')
    plt.show()

if __name__ == '__main__':

    main()
