# An intra farm stochastic model of Lumpy Skin Disease

This repository contains a stochastic, intra-farm model implemented using the package SimInf to investigate transmission dynamics of Lumpy Skin Disease. 

## Model formulation

We used a compartmental Ross-Macdonald model to reproduce the transmission dynamics between cattle and vectors (Figure 1A). The model considers a cattle population of 100 cows and a vector population of *Stomoxys calcitrans* with 50 flies per cow. 

The cattle population is subdivided into 10 compartments. The cattle population starts as fully susceptible ($S_C$). Upon infection, they begin their two-stage latent period (compartments $E_{C1}$ and $E_{C2}$). After the first stage of the latent period ($E_{C1}$, lasting $\frac{1}{a}$ days), infected cattle either follow a symptomatic route (with probability $p_S$), leading to the eventual expression of clinical signs, or an asymptomatic route (with probability $1-p_S$). Asymptomatic cattle complete their latent period in $E_{C2 Asymptomatic}$  after a duration of $\frac{1}{\gamma}$ days, before becoming infectious ($I_{C Asymptomatic}$) for a duration of $\frac{1}{\tau_{CA}}$  days and later recover ($R_C$). Future symptomatic cattle complete their latent period in $E_{C2 Symptomatic}$ after a duration of $\frac{1}{γ}$ days, before showing symptoms and becoming infectious ($I_{C Symptomatic}$). Symptomatic infected cattle can then either die from the infection (with a probability $p_D$) after a delay of $\frac{1}{α}$ days completed in $I_{C predeath}$, or recover with a probability ($1-p_D$) after a delay $\frac{1}{\tau_{CS}}$  days completed in $I_{C prerecovery}$. Due to the short time frame of the simulations (less than a year), recovered cattle remain immune to new infections in the R_C compartment, and non–disease-related demographic processes as well as animal movements are ignored. 

The vector population is subdivided into 2 compartments: susceptible $S_V$ and infectious $I_V$. Infectious vectors remain capable of mechanically transmitting the virus for a period of $\frac{1}{\tau_V}$  days following exposure, before becoming susceptible again. The population of vectors is considered to be at a stable equilibrium, with equal birth and death rates $\mu$.  
The force of infection exerted on a susceptible cow is given by

$$\lambda_C= \frac{k \cdot p_C \cdot  I_V}{N_C}$$

with $k$ being the biting rate of the vectors, $p_C$ being the probability that a cow becomes infected upon a bite by a contaminated vector, and $N_C$ being the total number of cows.
The force of infection exerted on a susceptible vector is given by 

$$\lambda_V=k (p_{VS}  \frac{I_{C Symp}+I_{C PreDeath}+I_{C PreRecovery}}{N_C} +p_{VA} \frac{I_{C Asymp}}{N_C} )$$

with $p_{VS}$ and $p_{VA}$ being the probabilities that a susceptible vector becomes contaminated if it bites a symptomatic infectious cow and an asymptomatic infectious cow, respectively.


## Interventions

