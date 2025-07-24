# Platt

See the [Platt website](https://cole-trapnell-lab.github.io/platt/) for a more comprehensive introduction. 


## Installation

Platt runs in the [R statistical computing environment](https://www.r-project.org/). It requires R >= 3.5.0. Platt is currently only available for Github install. 

### Required software

Platt builds on top of the [Hooke package](https://cole-trapnell-lab.github.io/hooke/install/). 

```r
devtools::install_github("cole-trapnell-lab/hooke")
```

Hooke builds on top of the [Monocle3 package](https://cole-trapnell-lab.github.io/monocle3/docs/installation/). 

```r
devtools::install_github("cole-trapnell-lab/monocle3", ref="develop")
```

Hooke depends on the [PLNmodels package](https://pln-team.github.io/PLNmodels/index.html).

```r 
remotes::install_github("pln-team/PLNmodels")
```

Finally, install the Platt package as follows: 

```r
devtools::install_github("cole-trapnell-lab/platt")
```
