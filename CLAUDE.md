# Teaching material bulk RNA-seq

## Goal

Write teaching material for students bioinformatics.
Use the style of `archive/old_examples/*.ipynb` (also use the same writing style)

The result should be a jupyter notebook.

## Task

Reproduce the results of Foxa1/Hnf4a overexpression in https://elifesciences.org/articles/73358
Focus only on analysing the bulk RNA-seq data presented in this publication.

The tutorial should at least cover all aspects in `archive/old_examples/*.ipynb`
Going from downloading the data -> mapping -> count matrix -> differential expresion
-> biological interpretation.

Data can be found here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182190

## Practical

If you need software, use mamba to install it. Keep environments minimal such
that dependencies do not clash

If you need a genome index, build it.

If something is unclear, or you are unsure about a decision, prompt the user.

Before starting write a plan and have it confirmed by the user