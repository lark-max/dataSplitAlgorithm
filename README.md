
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dataSplitAlgorithm

<!-- badges: start -->

This package can be used to split rainfall-runoff data into three data
subsets(Train, Test, Validation) with similar data distribution
characteristics. <!-- badges: end -->

## Installation

Here are two ways to install this package in Rstudio on your PC:  

Way 1:
Firstly ensure that you have installed the package `devtools`:  
```r
install.packages("devtools")
```
Then you can install the package like so: 
``` r
devtools::install_github("lark-max/dataSplitAlgorithm")
```
  
Way 2:
You can download the zip package provided in the `release` page,
and install the package by hand:  
```r
install.packages("dataSplitAlgorithm_1.0.0.tar.gz",repos = NULL,type = "source")
```
  
After installation,  
You can then load the package by:  
```r
library(dataSplitAlgorithm)
```

## Instruction

The package has several built-in data splitting algorithms, such as
SOMPLEX, MDUPLEX, SBSS-P, etc. These algorithms are encapsulated and
used by the user only by calling the function `dataSplit`.  
Specific call formats such as:

``` r
result = dataSplit(data,list(sel.alg = "SOMPLEX", writeFile = TRUE))
```

The parameters of the function are composed of two parts. The first
parameter `data` is usually the rainfall runoff data in data.frame or
matrix format, and set the first column as the subscript column. For the
specific format requirements, please refer to the documentation of the
built-in data set(`dataSplitAlgorithm_data_small`).

The second parameter `control` is a list of customized information, such
as the selected data splitting algorithm, the proportion of subsets,
whether the output file is required and the output file name, etc. These
information have default defaults, you can refer to the `par.default()`
provided in the package, which is shown as below:

``` r
[include.inp]   
Whether the input vector is included when calculating Euclidean distance between samples, which defaults to TRUE

[seed]  
Seed number, which defaults to 1000

[sel.alg]   
Algorithm selection, includes SOMPLEX, MDUPLEX, DUPLEX, SBSS.P, TC, which defaults to "MDUPLEX"

[prop.Tr]   
The proportion of data allocated to the training set, which defaults to 0.6

[prop.Ts]    
The proportion of data allocated to the test set, which defaults to 0.2

[Train] 
The name of the training set when the output file is required, which defaults to "Train.txt"

[Test]  
The name of the test set when the output file is required, which defaults to "Test.txt"

[Validation]    
The name of the validation set when the output file is required, which defaults to "Valid.txt"

[loc.calib] 
If using the TC algorithm, you can also specify the location of the splitting range, which defaults to c(0,0.6)

[writeFile] 
Whether to output the partition result to a txt file, which defaults to TRUE
```

## Example

The following example can be run directly from the user’s Rstudio
client.

``` r
library(dataSplitAlgorithm)
## basic example code
data("dataSplitAlgorithm_data_small")
result = dataSplit(dataSplitAlgorithm_data_small, control = list(sel.alg = "MDUPLEX",writeFile = FALSE))
#> [1] "Start the initial sampling..."
#> [1] "Initial sampling successfully!"
#> [1] "Start the loop sampling..."
#> [1] "Remaining unsampled data: 494"
#> [1] "Remaining unsampled data: 484"
#> [1] "Remaining unsampled data: 474"
#> [1] "Remaining unsampled data: 464"
#> [1] "Remaining unsampled data: 454"
#> [1] "Remaining unsampled data: 444"
#> [1] "Remaining unsampled data: 434"
#> [1] "Remaining unsampled data: 424"
#> [1] "Remaining unsampled data: 414"
#> [1] "Remaining unsampled data: 404"
#> [1] "Remaining unsampled data: 394"
#> [1] "Remaining unsampled data: 384"
#> [1] "Remaining unsampled data: 374"
#> [1] "Remaining unsampled data: 364"
#> [1] "Remaining unsampled data: 354"
#> [1] "Remaining unsampled data: 344"
#> [1] "Remaining unsampled data: 334"
#> [1] "Remaining unsampled data: 324"
#> [1] "Remaining unsampled data: 314"
#> [1] "Remaining unsampled data: 304"
#> [1] "Remaining unsampled data: 294"
#> [1] "Remaining unsampled data: 284"
#> [1] "Remaining unsampled data: 274"
#> [1] "Remaining unsampled data: 264"
#> [1] "Remaining unsampled data: 254"
#> [1] "Remaining unsampled data: 244"
#> [1] "Remaining unsampled data: 234"
#> [1] "Remaining unsampled data: 224"
#> [1] "Remaining unsampled data: 214"
#> [1] "Remaining unsampled data: 204"
#> [1] "Remaining unsampled data: 194"
#> [1] "Remaining unsampled data: 184"
#> [1] "Remaining unsampled data: 174"
#> [1] "Remaining unsampled data: 164"
#> [1] "Remaining unsampled data: 154"
#> [1] "Remaining unsampled data: 144"
#> [1] "Remaining unsampled data: 134"
#> [1] "Remaining unsampled data: 124"
#> [1] "Remaining unsampled data: 114"
#> [1] "Remaining unsampled data: 104"
#> [1] "Remaining unsampled data: 94"
#> [1] "Remaining unsampled data: 84"
#> [1] "Remaining unsampled data: 74"
#> [1] "Remaining unsampled data: 64"
#> [1] "Remaining unsampled data: 54"
#> [1] "Remaining unsampled data: 44"
#> [1] "Remaining unsampled data: 34"
#> [1] "Remaining unsampled data: 24"
#> [1] "Remaining unsampled data: 14"
#> [1] "Remaining unsampled data: 4"
#> [1] "MDUPLEX sampling complete!"
data("dataSplitAlgorithm_data_middle")
result = dataSplit(dataSplitAlgorithm_data_middle, control = list(sel.alg = "SBSS.P",writeFile = FALSE))
#> [1] "Total neuron: 60"
#> [1] "sampling on neuron: 1"
#> [1] "sampling on neuron: 2"
#> [1] "sampling on neuron: 3"
#> [1] "sampling on neuron: 4"
#> [1] "sampling on neuron: 5"
#> [1] "sampling on neuron: 6"
#> [1] "sampling on neuron: 7"
#> [1] "sampling on neuron: 8"
#> [1] "sampling on neuron: 9"
#> [1] "sampling on neuron: 10"
#> [1] "sampling on neuron: 11"
#> [1] "sampling on neuron: 12"
#> [1] "sampling on neuron: 13"
#> [1] "sampling on neuron: 14"
#> [1] "sampling on neuron: 15"
#> [1] "sampling on neuron: 16"
#> [1] "sampling on neuron: 17"
#> [1] "sampling on neuron: 18"
#> [1] "sampling on neuron: 19"
#> [1] "sampling on neuron: 20"
#> [1] "sampling on neuron: 21"
#> [1] "sampling on neuron: 22"
#> [1] "sampling on neuron: 23"
#> [1] "sampling on neuron: 24"
#> [1] "sampling on neuron: 25"
#> [1] "sampling on neuron: 26"
#> [1] "sampling on neuron: 27"
#> [1] "sampling on neuron: 28"
#> [1] "sampling on neuron: 29"
#> [1] "sampling on neuron: 30"
#> [1] "sampling on neuron: 31"
#> [1] "sampling on neuron: 32"
#> [1] "sampling on neuron: 33"
#> [1] "sampling on neuron: 34"
#> [1] "sampling on neuron: 35"
#> [1] "sampling on neuron: 36"
#> [1] "sampling on neuron: 37"
#> [1] "sampling on neuron: 38"
#> [1] "sampling on neuron: 39"
#> [1] "sampling on neuron: 40"
#> [1] "sampling on neuron: 41"
#> [1] "sampling on neuron: 42"
#> [1] "sampling on neuron: 43"
#> [1] "sampling on neuron: 44"
#> [1] "sampling on neuron: 45"
#> [1] "sampling on neuron: 46"
#> [1] "sampling on neuron: 47"
#> [1] "sampling on neuron: 48"
#> [1] "sampling on neuron: 49"
#> [1] "sampling on neuron: 50"
#> [1] "sampling on neuron: 51"
#> [1] "sampling on neuron: 52"
#> [1] "sampling on neuron: 53"
#> [1] "sampling on neuron: 54"
#> [1] "sampling on neuron: 55"
#> [1] "sampling on neuron: 56"
#> [1] "sampling on neuron: 57"
#> [1] "sampling on neuron: 58"
#> [1] "sampling on neuron: 59"
#> [1] "sampling on neuron: 60"
#> [1] "SBSS-P sampling complete!"
data("dataSplitAlgorithm_data_large")
result = dataSplit(dataSplitAlgorithm_data_large, control = list(sel.alg = "SOMPLEX",writeFile = FALSE))
#> [1] "Total neuron: 126"
#> [1] "sampling on neuron: 1"
#> [1] "sampling on neuron: 2"
#> [1] "sampling on neuron: 3"
#> [1] "sampling on neuron: 4"
#> [1] "sampling on neuron: 5"
#> [1] "sampling on neuron: 6"
#> [1] "sampling on neuron: 7"
#> [1] "sampling on neuron: 8"
#> [1] "sampling on neuron: 9"
#> [1] "sampling on neuron: 10"
#> [1] "sampling on neuron: 11"
#> [1] "sampling on neuron: 12"
#> [1] "sampling on neuron: 13"
#> [1] "sampling on neuron: 14"
#> [1] "sampling on neuron: 15"
#> [1] "sampling on neuron: 16"
#> [1] "sampling on neuron: 17"
#> [1] "sampling on neuron: 18"
#> [1] "sampling on neuron: 19"
#> [1] "sampling on neuron: 20"
#> [1] "sampling on neuron: 21"
#> [1] "sampling on neuron: 22"
#> [1] "sampling on neuron: 23"
#> [1] "sampling on neuron: 24"
#> [1] "sampling on neuron: 25"
#> [1] "sampling on neuron: 26"
#> [1] "sampling on neuron: 27"
#> [1] "sampling on neuron: 28"
#> [1] "sampling on neuron: 29"
#> [1] "sampling on neuron: 30"
#> [1] "sampling on neuron: 31"
#> [1] "sampling on neuron: 32"
#> [1] "sampling on neuron: 33"
#> [1] "sampling on neuron: 34"
#> [1] "sampling on neuron: 35"
#> [1] "sampling on neuron: 36"
#> [1] "sampling on neuron: 37"
#> [1] "sampling on neuron: 38"
#> [1] "sampling on neuron: 39"
#> [1] "sampling on neuron: 40"
#> [1] "sampling on neuron: 41"
#> [1] "sampling on neuron: 42"
#> [1] "sampling on neuron: 43"
#> [1] "sampling on neuron: 44"
#> [1] "sampling on neuron: 45"
#> [1] "sampling on neuron: 46"
#> [1] "sampling on neuron: 47"
#> [1] "sampling on neuron: 48"
#> [1] "sampling on neuron: 49"
#> [1] "sampling on neuron: 50"
#> [1] "sampling on neuron: 51"
#> [1] "sampling on neuron: 52"
#> [1] "sampling on neuron: 53"
#> [1] "sampling on neuron: 54"
#> [1] "sampling on neuron: 55"
#> [1] "sampling on neuron: 56"
#> [1] "sampling on neuron: 57"
#> [1] "sampling on neuron: 58"
#> [1] "sampling on neuron: 59"
#> [1] "sampling on neuron: 60"
#> [1] "sampling on neuron: 61"
#> [1] "sampling on neuron: 62"
#> [1] "sampling on neuron: 63"
#> [1] "sampling on neuron: 64"
#> [1] "sampling on neuron: 65"
#> [1] "sampling on neuron: 66"
#> [1] "sampling on neuron: 67"
#> [1] "sampling on neuron: 68"
#> [1] "sampling on neuron: 69"
#> [1] "sampling on neuron: 70"
#> [1] "sampling on neuron: 71"
#> [1] "sampling on neuron: 72"
#> [1] "sampling on neuron: 73"
#> [1] "sampling on neuron: 74"
#> [1] "sampling on neuron: 75"
#> [1] "sampling on neuron: 76"
#> [1] "sampling on neuron: 77"
#> [1] "sampling on neuron: 78"
#> [1] "sampling on neuron: 79"
#> [1] "sampling on neuron: 80"
#> [1] "sampling on neuron: 81"
#> [1] "sampling on neuron: 82"
#> [1] "sampling on neuron: 83"
#> [1] "sampling on neuron: 84"
#> [1] "sampling on neuron: 85"
#> [1] "sampling on neuron: 86"
#> [1] "sampling on neuron: 87"
#> [1] "sampling on neuron: 88"
#> [1] "sampling on neuron: 89"
#> [1] "sampling on neuron: 90"
#> [1] "sampling on neuron: 91"
#> [1] "sampling on neuron: 92"
#> [1] "sampling on neuron: 93"
#> [1] "sampling on neuron: 94"
#> [1] "sampling on neuron: 95"
#> [1] "sampling on neuron: 96"
#> [1] "sampling on neuron: 97"
#> [1] "sampling on neuron: 98"
#> [1] "sampling on neuron: 99"
#> [1] "sampling on neuron: 100"
#> [1] "sampling on neuron: 101"
#> [1] "sampling on neuron: 102"
#> [1] "sampling on neuron: 103"
#> [1] "sampling on neuron: 104"
#> [1] "sampling on neuron: 105"
#> [1] "sampling on neuron: 106"
#> [1] "sampling on neuron: 107"
#> [1] "sampling on neuron: 108"
#> [1] "sampling on neuron: 109"
#> [1] "sampling on neuron: 110"
#> [1] "sampling on neuron: 111"
#> [1] "sampling on neuron: 112"
#> [1] "sampling on neuron: 113"
#> [1] "sampling on neuron: 114"
#> [1] "sampling on neuron: 115"
#> [1] "sampling on neuron: 116"
#> [1] "sampling on neuron: 117"
#> [1] "sampling on neuron: 118"
#> [1] "sampling on neuron: 119"
#> [1] "sampling on neuron: 120"
#> [1] "sampling on neuron: 121"
#> [1] "sampling on neuron: 122"
#> [1] "sampling on neuron: 123"
#> [1] "sampling on neuron: 124"
#> [1] "sampling on neuron: 125"
#> [1] "sampling on neuron: 126"
#> [1] "SOMPLEX sampling complete!"
```

## Details

The package also integrates the function of adjunct validation. The
metric `AUC` is used to analyze the similarity of data distribution
features among the sample subsets obtained by various data segmentation
algorithms, which can be calculated by invoking the function `getAUC`.  
For example:

``` r
data("dataSplitAlgorithm_data_small")
res_split = dataSplit(dataSplitAlgorithm_data_small,list(sel.alg = "MDUPLEX",writeFile = FALSE))
res_auc = getAUC(res_split$Train,res_split$Validation)
```

In the above example, the value range of `res_auc` is \[0,1\], and the
closer the value is to 0.5, the more similar the data distribution
characteristics between the two sample subsets are.

## References

Chen, J., Zheng F., May R., Guo D., Gupta H., and Maier H. R.(2022).
Improved data splitting methods for data-driven hydrological model
development based on a large number of catchment samples, Journal of
Hydrology, 613.

Zheng, F., Chen J., Maier H. R., and Gupta H.(2022). Achieving Robust
and Transferable Performance for Conservation‐Based Models of Dynamical
Physical Systems, Water Resources Research, 58(5).

Zheng, F., Chen, J., Ma, Y., Chen Q., Maier H. R., and Gupta H.(2023). A
Robust Strategy to Account for Data Sampling Variability in the
Development of Hydrological Models, Water Resources Research, 59(3).
