# Bergman Kernels of Generalized Hartog's Triangles

### Overview
I have been informally working with Luke Edholm since my undergraduate at Michigan on some researh on 
several complex variables. This repo is highly undocumented as the research is on going and this exists
mostly as a sandbox to test various examples and notice patterns of coefficients in simpler cases.

The purpose of this repo was to teach myself some basic julia and use its computational power to generate 
examples and counterexamples as needed. Much of this code was duplicated from a previous repository written 
in MATLAB. This work is traditionally done in neither MATLAB nor Julia, but Sage or Mathematica. Unfortunately,
both have incredibly expensive licensing fees if you aren't tied to a university.


### The Research
Single-variable complex analysis is uniquely glorious in how well hard problems can become solvable. 
The Riemann mapping theorem is one of the great triumphs of complex analysis, providing
a way to classify all domains of the complex plane up to biholomorphism. Unfortunately, as you extend
this problem to the multivariate setting, the Riemann mapping theorem falls apart pretty immediately. However,
this task of classification can be done for specific families of higher dimension domains. We are interested 
in classifying an infinite family of so-called "Hartog's triangles" which are common places to find interesting
counter examples in SCV. We approach this problem by using a functional analysis tool called the Bergman Kernel 
which contains information about all of the holomorphic functions on a specific domain. Using functions to study
geometric or topological questions is common in complex analysis as holomorphic functions are particularly 
well-behaved. We have reduced this classification problem to proving certain polynomials do not vanish on the 
unit circle.

The PDF file left in the note gives a more specific overview about the state of the research and explicit 
formulas for the functions we are working with. Please note that the PDF is maintained elsewhere 
due to it being compiled in shared LaTeX files and therefore can be out of date. I try to update it whenever 
larger progress is made.

