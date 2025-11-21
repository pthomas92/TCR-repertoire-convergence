# TCR-repertoire-convergence
Basic repository to back up and distribute R script for calculating repertoire convergence stats (Jaccard Index) between immune states, and statistically enriched clones based on clonal over-representation (Fisher's exact test). Current focus is on differentiating abortive and lab-confirmed infection following SARS-CoV-2 human challenge, hence hard coded 'abortive' and 'symptomatic' variables in the data.

### Jaccard Index
Simple calculation of the jaccard index per donor and timepoint within each infection state, to allow for testing differential convergence between groups and timepoints. 2 abortive stats are calculated, one including transient infection and one excluding.

### Fisher's test
Since the over-representation test works at a single clonotype level, it usually must iterate through all clones and perform the statistical test, meaning compute times can be very long even with parallelisation. This script massively expedites this process by pre-computing the result based on possible sharing levels. E.g.,:
  1. There are 16 donors in the abortive infection group and 18 in the lab-confirmed group
  2. Pre-compute statistics based on all possible TCR sharing combinations. For example, pre-calculate the Fisher's exact test result for a TCR shared by 5/11 abortive donors and 0/18 lab-confirmed donors.
  3. Join these to the original data frame and identify enriched clones (p < 0.05).

P-values are adjusted at stage 3, since we have few donors and many TCRs, meaning post-hoc p-value adjustment inflates these values to 1.
