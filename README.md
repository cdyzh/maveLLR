# maveLLR
Calculate pathogenicity log likelihood ratios (LLRs) for MAVE datasets

## Installation
```R
devtools::install_github("jweile/maveLLR")
```
## Example usage

```R
#create some faux data representing our positive and negative reference variant scores
posScores <- rnorm(50,0.2,0.3)
negScores <- rnorm(30,0.8,0.2)

#calculate the LLR using kernel density estimation
llr <- buildLLR.kernel(posScores,negScores,bw=0.1,kernel="gaussian")

#draw a summary plot of the resulting LLR
drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)

#create some more faux MAVE scores
scores <- rnorm(500,0.5,0.3)
#transform them into their corresponding LLRs
result <- llr$llr(scores)
```
