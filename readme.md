Sample code for the Technical Lecture "DPD and Sparse Estimation" 
January 19th, 2022
Radio and Wireless Week 2022, Las Vegas, NV, USA.

Main scripts:
- FoM.m Figures of merit demonstrator
- dpd_DLA.m: DPD DLA demonstrator
- sparse_regression.m: Sparse regression demonstrator
    
Auxiliary scripts:
- Synthetic PA: 
    syntheticPA.m
    syntheticPA_inputoutput.mat

- Model generation:
    modeling.m: Modeling demonstrator
    buildX_GMP.m
    coeff_selection.m
    model_PA.m
    modelconfigGMP2.m
    omp_domp.m
    sel_indices.m
    
- Signal generation and analysis functions:
    signal_generation.m: 5G signal generator demonstrator
    generator5G.m
    analysis5G.m
    ACPR5G.m
    evm5G.m
    spectrum.m
    FFTinterpolate.m