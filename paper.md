---
title: 'PT_GSDesign: A SAS Macro for Group Sequential Designs with Time-to-event Data using the Concept of Proportional Time'
tags:
  - SAS
  - efficacy
  - error spending
  - futility
  - proportional time
  - sample size
authors:
  - name: Milind A. Phadnis^[Custom footnotes for e.g. denoting who the corresponding author is can be included like this.]
    orcid: 0000-0001-6472-9325
    affiliation: 1
  - name: Nadeesha Thewarapperuma
    orcid: 0000-0002-9306-2619
    affiliation: 1
affiliations:
 - name: University of Kansas Medical Center
   index: 1
date: 21 February 2021
bibliography: paper.bib
---

# Summary
Treatments that are found to be promising in a Phase II trial are studied more comprehensively in a Phase III trial where researchers 
aim to investigate the effectiveness and safety of the new treatment against the current standard-of-care. While traditional 
approaches require the calculation of a fixed sample size depending on the type I error, power and clinically important treatment 
effect, in the medical setting they suffer from the limitation that patients are continually being accrued into a 
study based on the accrual rate, the availability of qualified patients (based on inclusion/exclusion criteria) and the possibility 
of random dropouts among many factors. Thus, the primary outcome of interest is not available simultaneously on all patients 
and researchers may be interested to look at outcomes on the early enrollees and use that as a basis to decide whether 
the trial should be continued. Sequential testing in large-sized Phase III trials with interim points can be used to - {i} stop 
the trial early for overwhelming evidence of efficacy, {ii} stop the trial early for overwhelming evidence 
for futility, and {iii} continue the trial for lack of evidence of efficacy or futility.

A Group Sequential Design (GSD) formalizes the concept by providing a statistical framework under which either of the three decisions can 
be taken after looking at interim results. Ethical, financial and administrative requirements often guide the statistical 
designs of GSDs (see @enas; @jennison; @ellenberg). Such GSDs have been well developed for continuous and binary outcomes 
and have a long history starting with quality control applications (@wald) and progressing to the medical setting (@armitage). Vast literature 
is available on this topic in many books (@whitehead1; @jennison2; @pros; @dmi; @was) and overview articles 
(@whitehead2; @todd; @mazu). When dealing with time-to-event outcome, a repeated significance testing approach incorporating a 
family of designs (@pocock; @obrien; @wang) can be combined with the error spending method (@lan) to implement a GSD using a 
log-rank test or by using the proportional hazards (PH) assumption. Popular statistical software often implement GSDs for time-to-event outcome 
using the weighted and unweighted versions of the log-rank test either explicitly assuming exponentially distributed survival times or with the 
PH assumption and are able to incorporate complexities of survival outcomes such as random dropouts, prespecified accrual and follow-up times, 
varying accrual patterns, equal/unequal spaced interim testing points (looks), efficacy-only designs, efficacy and futility designs, 
binding and non-binding futility rules, and many other flexible features specific to time-to-event outcomes.
    
When the underlying assumptions that drive the analytical and simulation-based approaches using the framework of the log-rank test are not valid, 
hardly any alternate methods are available in literature or in standard statistical software. Recent developments in this field have considered 
relaxing the PH assumption in favor of a `proportionality of time (PT)' assumption leading to development of GSDs in the context of an accelerated 
failure time (AFT) model (@main). The authors have described various scenarios in the biomedical setting where their approach could be advantageous 
compared to the standard methods with the help of real-life examples. Their proposed GSD method provides an alternate approach when the PH assumption is 
not appropriate and allows various hazard shapes (increasing/decreasing monotonically over time, bathtub shaped, arc-shaped) using the generalized gamma 
distribution. The purpose of this paper is to present a fully functional SAS macro that can be used to implement their GSD method. 
The SAS macro incorporates multitude of design features specific to a two-arm GSD for time-to-event outcomes and includes validation for any parameters 
defined by the user, as well as suggestions for correcting erroneous input. 

# Acknowledgements

The High performance computing capabilities, which were used to conduct some of the analyses described in this paper, were supported in part by the National Cancer Institute (NCI) Cancer Center Support Grant P30 CA168524; the Kansas IDeA Network of Biomedical Research Excellence Bioinformatics Core, supported by the National Institute of General Medical Science award P20 GM103418; and the Kansas Institute for Precision Medicine COBRE, supported by the National Institute of General Medical Science award P20 GM130423.

# References
