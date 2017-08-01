[//]: # (Image References)
[ukf_ds1]: ./ukf_end_ds1.png
[ukf_ds2]: ./ukf_end_ds2.png
[nis_radar]: ./NIS_radar.png
[nis_lidar]: ./NIS_lidar.png

# Unscented Kalman Filter Project

This project implement an Unscented Kalman Filter to estimate the state of a moving object of interest with noisy lidar and radar measurements. 

![ukf_ds1]
![ukf_ds2]

# Current solution
The current solution have been developed using xcode 8.3.3 and tested with the term2 simulator provided, the reported RMSE is:   
Dataset 1: 0.06 0.08 0.30 0.37   
Dataset 2: 0.07 0.06 0.51 0.27   
The code is fairely organized, modular and it appears to run smoothly, basic optimization have been put in place.    

# Parameter tuning and NIS 

![nis_radar]
![nis_lidar]
