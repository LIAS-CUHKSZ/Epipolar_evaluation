Essential matrix estimation with the five point algorithm. Manolis Lourakis 2018

This code uses a five point algorithm in a RANSAC framework to compute an estimate of the essential matrix.
That estimate is subsequently refined by minimizing the symmetric distance from epipolar lines for RANSAC inliers.


Code uses the following functions by others (see links for origin):

- calibrated_fivepoint.m, calibrated_fivepoint_helper.c from http://vis.uky.edu/~stewe/FIVEPOINT (archived at https://web.archive.org/web/20170401223934/http://www.vis.uky.edu:80/~stewe/FIVEPOINT/)
  Note that calibrated_fivepoint_helper.c is a MEX file which should be compiled by typing (at the matlab prompt)
       mex calibrated_fivepoint_helper.c
  Alternatively, after editing essmatrix.m, you can also use
    -- calibrated_fivepoint_non_gb.m also from http://vis.uky.edu/~stewe/FIVEPOINT
    -- peig5pt.m from https://raw.githubusercontent.com/jianxiongxiao/SFMedu/master/peig5pt.m

- ransac.m from https://www.peterkovesi.com/matlabfns/Robust/ransac.m

- marquardt.m from http://www2.imm.dtu.dk/projects/hbn_software/marquardt.m (archived at http://archive.is/zEYmi)

