# ASNPnP_IRLS

The ASPnP algorithm enhanced with Iteratively Re-weighted Least Squares (IRLS).


## Running the code

To run the code in MATLAB, use the following command.

```
[R0, t0, w_final] = ASPnP_IRLS(U0, u0, K, w_prior)
```

Here, U0 is a 3-by-N matrix of 3D points, u0 is a 2-by-N matrix of their corresponding 2D projections, K is a 3-by-3 intrinsic camera matrix, and w_prior is an N-by-1 wector of prior weights to be assigned to each observation.


## Publications

Kindly cite the following paper if this code helps in your publications [paper-yet-to-appear].

Note: ASPnP was initially proposed [here](http://www2.maths.lth.se/vision/publdb/reports/pdf/zheng-kuang-etal-iiccvi-13.pdf). This code is only a slight, yet highly significant modification above ASPnP.
