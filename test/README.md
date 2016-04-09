##Test Results
- sampler_cal_product_test.py: This file is designated to test different implementation of the cal_product function in the sampler.py. The result shows the
original implementation by Shuo (loops) is significantly better than Chuqiao's numpy implementation. For simulation with real dimension (100 * 200000 * 20), Shuo's
method takes 130.41 sec while Chuqiao's method takes 193.15 sec.

- sampler_factor_test.py: This file is designated to test different numpy implementation when sampling the factor matrix. Mengqing implemented a method where numpy's dot product
is performed in every loop. Chuqiao's method, on the other hand, used extra space but extract the dot product out from the loop. Chuqiao's method is significantly better than Mengqing's method.

Last updated: by Chuqiao Ren on 04/08/2016 at 09:04 P.M.