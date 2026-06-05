Dear Dr Spurr,

I'm emailing with an update on your submission. Unfortunately, we are still waiting for the full set of review reports on your manuscript entitled "Environmental drivers shape the freshwater biofilm resistome across a national river network". In the meantime, we have decided to share the report(s) we have received. We hope you will find this useful but must stress that we have not made an editorial decision on your manuscript and ask that you do not provide any point-by-point response at this stage. Please see the reviewer report(s) below.

Thank you for being patient with us.

Best regards,

Danielle Troppens
Editor
Nature Communications



Reviewer’s Comments:



Reviewer #1 (Remarks to the Author)

The study “Environmental drivers shape the freshwater biofilm resistome across a national river network” by Spurr et al. evaluates the resistome abundance and identifies distinct resistome signatures in freshwater biofilm across riverine network nationwide under anthropogenic and hydro-climate events. The study is intensive and well put together and contains exciting insights about the resistome abundance and these findings directly support the development of an evidence base for future regulatory processes and move towards a risk-based approach to managing AMR in the environment. It further advances the prospects for environmental AMR surveillance by highlighting the relative contributions to the environmental resistome from multiple stressors.
"The core concepts of this study are quite promising and address an important topic. However, to fully realize the potential of the work, the manuscript requires further refinements. I believe a major revision would significantly strengthen the paper for future consideration."
Major Comments:

1. Authors have not clearly mentioned how the microbial community structure of the ARGS containing host is changing across the sewage, agriculture, and rainfall clusters.
2. A network-based approach provides a robust framework for evaluating the distribution of ARGs, allowing for a deeper understanding of how riverine connectivity influences their associations and dissemination.
3. MAGs based recovery of multiple uncultured ARG-hosting bacteria will further highlights the limitations of traditional culture-based surveillance and the importance of including yet-uncultivated microbes in AMR monitoring frameworks. These taxa may harbour resistance genes that have not yet been described and may act as vectors for ARG transfer.
4. Authors have performed the ARG analysis with contigs which will not shed light on ARGs residing in the microbial dark matter and raise concerns about their role as previously overlooked resistance reservoirs.
5. The study would benefit from correlating ARG abundance with quantified antibiotic levels within the riverine biofilms as well as all specific stressors like sewage, agriculture, etc..
6. Authors have mentioned in the Introduction part that “Rivers, acting as drainage basins for terrestrial and urban landscapes, serve as primary receptors for a complex mixture of antimicrobials, heavy metals, and resistant bacteria that originate from anthropogenic activities”, therefore, to understand and find an effective strategy to manage antibiotic resistance, it is necessary to analyze ARGs along with metal resistance genes and mobile genetic elements because ARGs and metal resistance genes co-occur due to shared mode of action, which could be co-selected and spread together.

Reviewer #2 (Remarks to the Author)

The paper by Spurr and colleagues presents a national survey of antimicrobial resistance determinants (AMDs) in metagenomes derived from biofilm samples from rivers across England, with repeat samples taken across different seasons. The study is large in scale and potentially impactful. The focus on biofilm samples rather than water samples makes sense and the conclusion that models of AMD and pathogen dynamics in waterways should include hydrologic data are well justified. However, we suggest several improvements that would strengthen the paper.

Major comments
1. The major interpretation that rainfall “mobilizes” AMDs is oversimplified. The study bases this statement on the abundance of AMDs per genome, which is ultimately a measure of an increase in relative abundance of taxa carrying a higher number of AMDs per genome in response to rainfall. However, the relative abundance of these genes/organisms can’t just be understood as a function of “mobilization” without watershed-level studies and consideration of absolute abundance (per mass or volume). An increase in relative abundance of AMDs per genome as measured here is a result of several factors including mobilization, colonization, and/or selection and this cannot be compared/contrasted to the dilution hypothesis, which focuses on the concentration of pollutants/antimicrobials in water. In our view, mobilization could be directly measured as a time course within a watershed where there would be an increase in prevalence and/or absolute abundance of AMD-carrying microbes at downstream sites in response to rainfall. The paper has a little more enlightened discussion of factors that could lead to increased abundance of AMD-carrying microbes in the discussion (lines 576) because this section mentions runoff of nutrients or mobilization of organisms, but mobilization of antimicrobials followed by selection for AMDs should also be discussed. In this vein, the paper discusses the “pulse dynamic” model, which involves scouring surfaces and resuspending cells into the water column, but if surfaces are scoured then that would lead to a decrease in total abundance in biofilms. This is just another example of where the logic needs to be tightened.
2. The methods are inadequate particularly regarding sampling and statistics. The paper cites another paper for information on sampling methods but the paper is very focused on sampling a different environment (rather than water) and this should be better described, particularly how the biofilm samples were collected. More importantly, the methods mention several statistical tests (some parametric and some non-parametric) but the paper almost never describes which analysis/figure uses which statistical test and so the reader/reviewer must trust that statistical tests are justified. In this regard, the figure legends are generally inadequate. Additionally, when describing bioinformatic tools and software, all parameters and settings used must be specified.
3. In general, it’s a little strange to base the whole study on contigs and not MAGs, especially since other papers from this research group analyzed MAGs from these same samples. We think this is possibly due to the relatively low number of MAGs possible due to a combination of (i) complex community structure, (ii) the shallow nature of the metagenomes, and (iii) the difficulty binning AMDs, which are often on mobile genetic elements. But the mobility of AMDs could also complicate assignment of contigs to specific taxa, which is not discussed. We personally would prefer an analysis of AMDs in MAGs and on contigs but we don’t think it’s absolutely necessary. But minimally the paper should deal with these issues in some way.
4. Based on our experience and recent publications (Comprehensive taxonomic identification of microbial species in metagenomic data using SingleM and Sandpiper | Nature Biotechnology) we’re quite skeptical of the accuracy of Kraken2, especially because the methods don’t describe a contig size cutoff and because of the prevalence of many AMDs on MGEs, which are probably hard to assign to specific organisms. We suggest a couple of ideas that could increase confidence in these assignments: (i) use of another approach such as SingleM or (ii) an analysis requiring contigs of certain sizes (e.g., 10 kbp). If these give similar results this would increase confidence in the assignments.
5. The paper doesn’t describe filters for mitochondria or chloroplasts, only for reads mapping to the human genome. Are contigs from mitochondria and/or chloroplasts being included in the Pseudomonadota and Cyanobacteriota counts and if so, how is this affecting the analysis? These two phyla are the most abundant in the results, and there is a large amount of discussion about both.
6. The paper requires a better treatment of gene acronyms. They should be clarified at first use and there should be a list of all 114 gene acronyms with descriptions, ideally with an accounting of each gene in each sample and its assigned host(s). This could be a supplemental table. In general, a lot of data would be needed to provide the support and transparency demanded by NPG.
7. Writing needs quite a lot of work. Below are some themes we observed, but in general one of the authors with a strong command of writing and grammar should go through the paper carefully. (i) inconsistent use of Oxford commas. (ii) consistent references to specific taxonomic names while clarifying the rank after the taxon name (e.g., the Bacteroidota phyla). These are restrictive appositions where the more general rank needs to come before the more specific taxonomic name and the number needs to agree (e.g., the phylum Bacteroidota). Consider “my Susan friends” versus “my friend Susan.” (iii) many sentence fragments beginning with subordinate conjunctions (e.g., whereas, whilst, etc.).

Minor comments
1. Line 41: is this death toll from HIV AND malaria combined or separate?
2. Which alignment tool did you use within RGI? DIAMOND or blast? Please specify and cite accordingly.
3. Line 133: what are “Kraken2 taxonomies”? Kraken2 uses NCBI taxonomy by default. Please specify (and correct pluralization if referring to a single taxonomy).
4. Line 145: here the logic behind “ARG per taxon” is described, where the abundance of ARGs was normalized to the genus level. We don’t understand why the phrase “ARG per genus” isn’t used instead, as that would be more accurate and specific.
5. Lines 198-199: What is evidence of this? Distinct from what? Would need support to say this.
6. “We observed a robust positive linear correlation between total contig coverage and ARG abundance that did not plateau (Figure 1b), suggesting that AMR load may be underestimated by the sequencing depth used in this study.” We agree, but Figure 1b should be improved by actually demonstrating this linear relationship statistically rather than relying on observation.
7. What does this sentence mean: “This also suggested that the sequenced resistome scaled density-dependently and was consistent with the overall microbial biomass of these biofilm communities.”
8. The histograms in figure 1a and 1c have unclear bin sizes. Please revise these panels.
9. Lines 229-230: what does “sorted from most-least Pseudomonadota and least-most Cyanobacteriota” mean?
10. Figure 2 is good but none of the results are significant so maybe this would be better placed in the supplement.
11. Figure 3 and onward (presumably) use Pearson correlation coefficient. Please clearly justify the use of this test. The normal distribution is shown in Figure 1c, but the linear relationship shown in Figure 1d needs to be demonstrated statistically.
12. Line 301: But why is the correlation between diaminopyrimidine resistance and woodland cover in the South West and central Northern England interesting? Elaborate or rephrase.
13. Line 306: “This meant that the driving environmental factors could not be easily identified as many machine learning approaches did not model well or yield clearly defined clusters.” What does this sentence refer to? The paper doesn’t describe any machine learning. This is unacceptable.
14. Line 343: “strong, yet distinct” these two adjectives are in agreement with each other, not conflict.
15. Line 428: remove “abundances” (and pluralize ARGs).
16. Line 432: “(formerly Proteobacteria)” should be removed.
17. Line 437: What in particular is so interesting about these ARGs being hosted by the phylum Nitrospirota? Elaborate or rephrase.
18. What is “Genera hosting ARGs within phyla”? This makes no sense. All genera belong to a phylum, right?
19. The paper often discusses ARGs on various taxa, which doesn’t make sense. Then in some cases they are discussed in various taxa, which sounds better. In the sentence below, the terminology is switched. This should be clarified. “in the Pseudomonadota phylum, primarily on Erwinia spp., but also on…”
20. Line 486, Flavobacterium is not a cyanobacterium. Also, this genus is not mentioned in the section so why should it be mentioned in the conclusion sentence for the section?
21. Lines 498 & 502: should be “Enterobacter spp.” not sp.
22. Line 506: where does “15,737 species-ARG detections” come from? This is the first anything like this is mentioned. The methods and main text should be improved to make the differences between analyses using “ARGs per taxon [genus]” and analyses using these “species-ARG detections” more clear to the reader. Supplemental data table(s) like suggested above could also make this clearer.
23. Line 523 refers to other large-scale studies similar to this. It would be better to mention some.
24. Lines 572-573: The correct name of the phylum is “Cyanobacteriota,” not “Cyanobacteria.”
25. Some figures have problems with the text bleeding into the figure and should be corrected, and the text-heavy figures should be generally improved for legibility.
26. Data Availability: the ENA accession number for this data should still be included in this section, even if it was previously published.
27. Data Availability: there is a Zenodo link to the “processed AMR data linked by the sample metagenome ID,” but it is currently embargoed so we are not able to review it.
28. Code Availability: section is missing—presumably there is a Github page with the Snakemake workflows and other code described in the methods?

