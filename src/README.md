#**Extended Kalman Filter**


**Extended Kalman Filter Project**

[//]: # (Image References)

[image0]: ./../results/best_result_dataset1.png "best_result_dataset1.png"
[image1]: ./../results/best_result_dataset2.png "best_result_dataset2.png"



Of course the Kalman filter based on the sensor fusion of Lidar and Radar produces the best results.
I tried out using just Lidar and Radar (actually this just means commenting out line 194 or 199 in file FusionEKF.cpp) and both results are worse compared to the sensor fusion Kalman filter. 
The simple reason is that more measure points are feed into the Kalman filter compared to using only one sensor signal.
The other reason is that Lidar and Radar can balance their individual strenghts and especially weaknesses.

####Advantages of Lidar:
* spacial resolution

####Disadvantages of Lidar:
* cannot measure velocity directly
* depending on weather conditions


####Advantages of Radar:
* high distance
* independent from weather conditions

####Disadvantages:
* Low resolution (especially in vertical direction)
* Reflection on static objects

Setting the noise_ax and noise_ay both to value 50 produces the best results - see the following two screenshots

![best_result_dataset1.png][image0]

![best_result_dataset2.png][image1]