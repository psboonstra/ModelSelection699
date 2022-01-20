# ModelSelection699
Materials for 699 Lecture on Model Selection

Go to www.umich.edu/~philb/Winter2022Slides/ModelSelection.html to view the slides

- `ModelSelection.R` is all of the R code I used to create the presentation. 
It is the result of running `knitr::purl("ModelSelection.Rmd")`, where 
`ModelSelection.Rmd` is a markdown file (that I have not shared) that I wrote

- `stepAICc.R` is a simple tweak to the `stepAIC()` function that implements
the corrected AIC. It contains two functions: `extractAICc()` and `stepAICc()`. 
Currently this only works for logistic regression; the proper approach here would
be to create a generic `extractAICc()` function just like `extractAIC()`

- `varselection_sim.R` does the simulation study from the beginning of the lecture

- `bdiag.csv` is the Breast Diagnosis data from Street, Wolberg, and Mangasarian (1993). 
See the presentation for the full citation
