# Installation

Platt runs in the [R statistical computing environment](https://www.r-project.org/). It requires R >= 3.5.0. Platt is currently only available for Github install. 

> **_NOTE:_** Platt is currently in the beta phase of its development. The documentation on this page is also still under construction. Not all features currently implemented have been completely documented. Please report any issues to your [github page](https://github.com/cole-trapnell-lab/platt/issues). 


### Required software

Platt builds on top of the [Hooke package](https://cole-trapnell-lab.github.io/hooke/install/). 

```r
devtools::install_github("cole-trapnell-lab/hooke", ref="develop")
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

```r
devtools::install_github("cole-trapnell-lab/platt", ref="develop")
```
