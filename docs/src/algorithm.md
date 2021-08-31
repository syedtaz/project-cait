# Algorithm

Currently the package supports a delayed stochastic simulation system using a modified next reaction method to model the zebrafish MRNA concentrations. Please see Ay & Ozubdak[^1] for more details on the research project -- this page only contains a brief overview.

## Next Reaction Method

Please see Anderson[^2] for a theoretical and language agnostic discussion on the algorithm implemented in this package. 

### Terminology

- ``t`` refers to the current time in the system
- For any reaction ``k``
    - ``P_k`` refers to the propensity of the reaction
    - ``T_k`` refers to the internal time  the reaction
    - ``a_k`` refers to the propensity function  the reaction
    - ``r_k`` refers to the ``~\mathcal{U}(0,1)`` random number generated for  the reaction
    - ``S_k`` refers to the completion time channel of the reaction (only if it is a delayed reaction)
- ND refers to reactions that have no delays. Their effects are applied immediately.
- CD refer to reactions that have a delay but their change is applied when the reaction is completed.
- ICD refers to reactions that change the state of the system both when they initiate and when they complete.

### Pseudocode

1. Set the initial number of reactants and set ``t=0``. For every reaction ``k``, set ``P_k`` and ``T_k`` to zero and initialize ``S_k = [\infty]`` if it is a delayed reaction.
2. For each reaction ``k``,
    1. Calculate the propensity functions ``a_k``.
    2. Generate ``r_k`` and set ``P_k = \ln{\dfrac{1}{r_k}}``.
    3. Set ``\Delta t_k = (P_k - T_k)/a_k``
3. Set ``\Delta = \min{(\Delta t_k, s_k(1) - t)}``
4. Increment ``t`` by ``\Delta``.
5. Let ``\mu`` be the reaction we chose. 
    - If we are completing a ICD reaction apply the effect of the reaction and delete the first element of ``s_{\mu}``.
    - If we are initiating a non-delayed reaction, apply the effect of the reaction immediately.
    - If we are initiating a CD reaction, update ``s_{\mu}`` by inserting ``t + \tau_{\mu}`` into ``s_{\mu}`` in the second to last position.
    - If we are initiating an ICD reaction, apply the initiation effects of the reaction and update ``s_{\mu}`` by inserting ``t + \tau_{\mu}`` into ``s_{\mu}`` in the second to last position.
6. For each reaction ``k``, set ``T_k = T_k + a_k \cdot \Delta``. 
7. If reaction ``\mu`` was initiated, let ``r`` be ``~\mathcal{U}(0,1)`` and increment ``P_{\mu}`` by ``\ln{\dfrac{1}{r}}``.
8. Go to step 2, or quit.

## Rejection / Direct Method

These algorithms were tried out but were not generally as fast or readable as NRM. Stable and efficient versions might be added in the future.

[^1]: Ay A, Knierer S, Sperlea A, Holland J, Özbudak EM. Short-lived Her proteins drive robust synchronized oscillations in the zebrafish segmentation clock. Development. 2013 Aug;140(15):3244-53. doi: 10.1242/dev.093278. PMID: 23861061.
[^2]: J. Chem. Phys. 127, 214107 (2007); https://doi.org/10.1063/1.2799998