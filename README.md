# COXEN (CO-eXpression ExtrapolatioN)
An implementation of and analysis with the COXEN algorithm.

## Introduction & Purpose
This was originally implemented in our lab at the Animal Cancer Center (Colorado State University) in R by Jared Fowles based on [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2831138/), which predicts patient survival based on microarray gene expression. I expanded on the R implementation, automating certain parts, applying and testing additional statistical models (including some machine learning models) and creating the "reverse" implmentation (choosing a chemotherapy based on gene expression). 


See our publications for more information:
  * [Intra- and interspecies gene expression models for predicting drug response in canine osteosarcoma](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4759767/) 
  * [Drug Selection in the Genomic Age: Application of the Coexpression Extrapolation Principle for Drug Repositioning in Cancer Therapy](https://www.liebertpub.com/doi/10.1089/adt.2015.29012.dlgdrrr)


This repository contains an implementation of COXEN in R with some of the original scripts. The purpose of this repository is to update the implementation for usability as this was mainly used in just our lab and I am unable to find the original implementation online. 

## Usage
These are currently for use in our lab and as such may not work as well for others. Once the implementation is cleaned up and more generalized I will add a description of how to install and use this code. For now you may explore the code and run using R. 

## Plans
The remaining work to be done on this code is to improve the efficiency and clean the code up. This code was written in 2013-2014 when I first started learning R and I've learned a lot from originally writing it and since then on other projects. Once it is a bit more cleaned up, perhaps released as a package in R, I will also attempt to implement it in Python for practice.

## Contributing
If you would like to contribute, feel free to create an issue or open a pull request. 

## Licensing
Jae K. Lee and Dan Theodorescu have intellectual property rights in the COXEN principle. We merely coded the vague idea with our own interpretation for research purposes. Thus, if this code is used, modified, repurposed, etc for research it should reference the original [publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2831138/) and the papers listed above since this code was created for that work. I have added an MIT license for the code, but if this is incorrect because of these details, I would love to know to update it. 
