---
title: 'substrateminer: A Python package to investigate protein substrate repertoires'
tags:
  - Python
  - visualisation
  - enzyme
  - proteolysis
  - substrates
  - bioinformatics
authors:
  - name: Robert QIAO
    orcid: 0000-0003-4569-6921
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
  - name: School of Biological Sciences, Flinders University, Bedford Park, SA 5042, Australia
    index: 1
  - name: Digital Research Services, Flinders University, Bedford Park, SA 5042, Australia
    index: 2
date: 29 Aug 2024
bibliography: paper.bib
---

# Abstract
`Substrateminer` is a Python package that provides a suite of investigative tools to develop protein substrate repertoires for a given protease of interest.  `Substrateminer` provides tools for analysing protease cleavage consensus and visualising the conservation and strength of consensus. `Substrateminer` is a complete suite that enables potential substrate identification for a protease of interest in a target organism and extended functions allow for further discovery of associated biological and pathological functions of detected substrates.

# Introduction
Enzyme hydrolysis is a fundamental biochemical process in biology, and underpins many biological processes in nature from food digestion in humans [@campbell2005biology] to renewable energy production with cellulosic ethanol [@GUO2023128252;@osti_6943966]. Proteases are a class of enzymes that carry out such hydrolytic functions in biology by cleaving peptide bonds. In human physiology, a wide variety of critically important biophysiological processes depend on well-regulated protease functions, including protein degradation, cell signalling, and immune response. [@R800035200;@Bond2019;@ajplung001562016] Consequently, the inhibition of target proteases has long been one intervention strategy applied clinically to treat diverse pathophysiological conditions including viral infections like HIV, hepatitis-C; metabolic dysfunctions like type 2 diabetes [@Scott2017]; and increasingly, cancers [@Manasanch2017]. In the wake of the COVID-19 pandemic and increased respiratory diseases, proteases have emerged as a major therapeutic target for curbing fatality due to viral infections.[@ijms25158105;@v16030366;@ijms22115762]

Nevertheless, the protease network in many organisms especially mammals including humans is complex and heavily-entangled, and the regulatory scope for many target proteases is still not well-understood. The ability to gather and identify novel protease substrates is an important step towards a better understanding of the biological functions of proteases and accelerate future developments of new protease inhibitors for downstream therapeutic applications. Here, the introduction of `substrateminer` aims to provide a comprehensive suite of tools to streamline the investigation and develop protein substrate repertoires for a given protease of interest based on a consensus algorithm.

# Methods
`Substrateminer` consists of three main categories of functions (\autoref{fig:workflow}), namely consensus investigation (i.e., with modules msa and consensus), substrate mining (i.e., with module miner) and substrate biopathological pathway searching (i.e., with module pathfinder).

![Substrateminer workflow\label{fig:workflow}](graphics/workflow.png)


A typical workflow for adopting `substrateminer` to investigate protein substrate repertoires for a given protease often follows three main steps:

## Consensus investigation
For a target protease with a few known cleaving substrates, finding a cleavage consensus is often beneficial. `Substrateminer` provides access to a mechanism to visualise the cleavage consensus motif and calculate the strength of conservation at each site along the motif. The strength of conservation is calculated on entropy and `substrateminer` is also packaged with auxiliary tools to provide multi-sequence alignment (MSA) access to help prepare the input data for consensus investigation. Consensus visualisation is based on the implementation of weblogo [@Crooks01062004]. (\autoref{fig:weblogo})

![Consensus investigation - Visulisation Example\label{fig:weblogo}](graphics/weblogo-2.png)

## Substrate mining
Upon the identification of a cleavage consensus, `substrateminer` provides a mechanism to mine potential substrates based on three filter strategies:

1. _Cellular location_: filter potential substrates based on the cellular localisation of the intended target. For example, if the target protease is a secreted protease, the substrates mined should be extracellular proteins.
2. _Target Size_: filter potential substrates based on the size of the intended target. For example, if the target substrate is a small protease, the substrates mined should be small proteins or peptides to a size of determined value.
3. _Consensus_: filter potential substrates based on the conservation of the cleavage consensus. At the current release, both endopeptidase and exopeptidase cleavage consensus are supported.

## Substrate pathobiological pathway searching
To further investigate the pathobiological relevance of novel substrates, `substrateminer` implemented a search strategy to retrieve known biological processes and disease associations based on the _de facto_ standard KEGG database [@10.1093;@pro3715]. In the present release, the `pathfinder` module in the `substrateminer` provides three main functions namely KEGG accession conversion, biological processes tracing and human disease mapping, where KEGG accession conversion enalbes a gateway to KEGG databases where further information including enzyme nomenclature and diseases-related network databases can be easily accessible. A more comprehensive implementation of KEGG access is scheduled for future releases.

# References