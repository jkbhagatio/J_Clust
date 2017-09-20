# J_Clust
A semi-automated, MATLAB-based spike-sorting software package for tetrode recordings of extracellular brain regions. Outperforms similar packages including M-Clust, WaveClus, XClust and SNAP Sorter.

![j_clust2](https://user-images.githubusercontent.com/14895866/30189822-6bc536c6-9404-11e7-8c08-02670035f634.jpg)

Manual: [J_Clust2_Manual.pdf](https://github.com/jaib1/J_Clust/files/1297605/J_Clust2_Manual.pdf)

Please see the Wiki in this repository for introductory videos on using J_Clust2, FAQs, planned future development and helpful references.

If you have any questions, comments, or concerns, or would like any additional information, please contact me - Jai Bhagat - at jaib1@mit.edu

This code is distributed and protected under the Fair Source License. For more information, see [here](https://fair.io/)

**Note on computation: Currently, on a system that runs 16 GB RAM, if loading in raw signal, J_Clust is optimized to handle < 100000 spikes per session, and OPTICS is optimized to cluster < 40000 spikes. If loading in pre-detected spike waveforms, or for systems with > 16 GB RAM, the number of spikes the program can handle is larger.
