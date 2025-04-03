# Code and data associated with: A survey of mammal and fish genetic diversity across the global protected area network
### Schmidt et al. Conservation Letters (2025) 
https://conbio.onlinelibrary.wiley.com/doi/full/10.1111/conl.13092

This folder contains:
- Scripts for data preparation and analysis are stored in order (1-4) for terrestrial and marine data.
- Data folder: datasets used for analysis.
- Networks folder: networks used for spatial regressions.

Data description: <br>
<b>genetdata_PA_terr</b> is the dataset used for terrestrial analyses. <br>
<b>genetdata_PA_mar</b> is the dataset used for marine analyses. See also Karachaliou et al. https://doi.org/10.1111/mec.17711 <br>
<br>
Column descriptions:
- filename: name of associated genotype file
- pop: unique sample site ID
- species: species binomial name (Genus_species)
- author: author of original dataset
- genus: species genus
- data_doi: DOI of original dataset
- paper_doi: DOI of original publication
- num_individuals: number of individuals at site
- num_loci: number of individuals at site
- gene_diversity: gene diversity at site
- allelic_richness: allelic richness at site rarefied to minimum 10 alleles
- allelic_richness_rbysp: allelic richness at site rarefied to minimum number of alleles sampled per species
- global_fst: population-specific FST (Weir & Goudet 2017 Genetics)
- Ne: contemporary effective population size
- Ne_lower & Ne_upper: 95% confidence intervals of Ne estimate
- PA_bin: protected status of site; 0 = not in protected area, 1 = within protected area
- PA_dist_neg: distance to nearest edge of a protected area. Value is negative when a site is located inside a protected area
- PA_count_10km, _50km, _100km, _150km: number of protected areas within 10, 50, 100, and 150km buffers around each site
- n_pops: number of sites per dataset
- PA_prop: proportion of sites located inside protected areas per dataset
- WDPA_PI: ID of protected area
- AREA_KM: area of protected area
- year: year of enactment of protected area status
- IUCN_CA_num: IUCN protected area category (integers 1-6; 1 = restrictive)
- adult_mass_g: average adult mass of species in grams
- lon, lat: site longitude and latitude. Coordinates have been removed for sensitive species and datasets. All coordinates have been truncated to 1 decimal digit to protect location data but enable analysis.
- total_alleles: total allele count per dataset
- pa_allele_count: count of alleles sampled within protected area per dataset
- prop_alleles_inPA: proportion of alleles sampled within protected area per dataset
