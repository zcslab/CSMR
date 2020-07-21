## Supervised Learning within Different Subspaces using Mixture Modeling

### Motivation

A tremendous amount of applications deal with relating a response variable to a set of covariates. The homogeneity assumption that the relation maintains for the entire population is often inadequate, when there exists complex subspace architecture in the data. Detection of all subspace domains comes with inhibitive computational cost. In this work, we proposed an efficient supervised subspace learning framework, which are particularly useful in detecting the dynamic dependencies between subsets of features and a response variable specific to unknown subspaces in the input data, namely **C**omponent **S**parse **M**ixture **R**egression (**CSMR**). We demonstrated CSMR to be very promising to deal with high dimensional data with complicated unknown subspace structures. Experimental evaluation on simulated benchmark data demonstrated that CSMR can accurately identify the subspaces on which subsets of features are explainable to the external variable. Application of CSMR on one real-world cancer cell line drug sensitivity data demonstrated that the subsets of feature and subspaces of subject identified by the proposed method have more favourable performances.

![image](https://github.com/zcslab/CSMR/blob/master/img/CSMR_frame.png)

## Installation

```
install.packages("devtools")
devtools::install_github("zcslab/CSMR")
```

## Usage
```
result=CSMR(x,ync,max_iter)
```

## Arguments
* ```x``` High dimensional data matrix. The row is the object and the column is the feature.
* ```y``` The external supervised variable.
* ``` nc ``` The component number in the mixture model.
* ``` max_iter ``` The maximum iteration number.


## Value
Result list contains five elements: ```coffs``` shows the coefficient matrix in mixture regression model; ```clus``` is the predicted cluster membership for each object; ```x``` is the input high dimensional data matrix; ```y``` is the input external supervised variable; ```yhat``` is the predicted external variable based on the mixture model.

## Example
```
library(CSMR)
X = example_data$x
y = example_data$y
result=CSMR(X,y,nc=2,max_iter=50)
```

### Questions & Problems

If you have any questions or problems, please feel free to open a new issue [here](https://github.com/zcslab/CSMR/issues). We will fix the new issue ASAP.  You can also email the maintainers and authors below.

- [Wennan Chang](https://changwn.github.io/)
(wnchang@iu.edu)

PhD candidate at [Biomedical Data Research Lab (BDRL)](https://zcslab.github.io/) , Indiana University School of Medicine
