# MSC

The Multiscale Coupled Model - Matlab Implementation  

For questions or bugs please open an issue or report to sebastiano.stipa@vki.ac.be

This model has been implemented by Sebastiano Stipa, during his PhD at the University of British Columbia. The MSC model couples a wake model to 
an atmospheric perturbation model in order to capture the wind farm interaction with the atmospheric boundary layer (ABL). The MSC model has been 
published in the Journal of Wind Energy Science (https://wes.copernicus.org/articles/9/1123/2024/wes-9-1123-2024.html) where its theory is explained 
and verification is carried out against large-eddy simulations conducted using the Toolbox fOr Stratified Convective Atmospheres 
(TOSCA, https://github.com/sebastipa/TOSCA).  

This implementation also feature a deep-array module which, upon activation, allows to take into account the evolution of cluster wakes better than 
wake models alone. This part of the model has not been published in any journal yet, but it is described in Sebastiano Stipa's PhD thesis (https://open.library.ubc.ca/soa/cIRcle/collections/ubctheses/24/items/1.0442080). 

The MSC model couples the Bastankhah and Porté-Agel (2014) wake model (https://www.sciencedirect.com/science/article/abs/pii/S0960148114000317) to Allaerts and Meyers' (2019) three-layer model 
(3LM, https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/sensitivity-and-feedback-of-windfarminduced-gravity-waves/A8C7C0206202618C01FD6FBFC2BEC87F), with 
optional deep-array effects modeled using a largely extended version of Meneveau's (2012) top-down model (TDM, https://www.tandfonline.com/doi/full/10.1080/14685248.2012.663092). 

The innovation of the MSC model lies in the coupling strategy between these three models, which operate at different spatial scales. In particular, the MSC coupling strategy allows to eliminate double counting effects between the different models, thereby enabling a transfer of information between the different scales that allows to take into account the most 
important aspect of the flow around wind farms, namely individual wakes, local blockage effects, global blockage effects and deep array effects. 

A special thanks in developing the model goes to Joshua Brinkerhoff, Nick Robinson, Arjun Ajay and Dries Allaerts who all added, especially Dries, valuable intellectual contributions 
to the ideation and development of the MSC model. 

The MSC model has been already incorporated into commercial wind farm design codes (OpenWind), as well as comprehensive open-source packages (WAYVE, https://iopscience.iop.org/article/10.1088/1742-6596/2767/9/092079). This implementation is the original MSC implementation, which stands here as a proof of concept for the model and as a platform for future 
implementations and additions. 

We hope that you will find the model useful. 

## Running the model 

The MSC model is fairly easy to run, in that only a working Matlab installation is required. In order to run the model, copy the content of the ``src`` directory into the location 
where you wish to run the model, so that the repository folder is left *as-is*. The name of the new location that is a copy of the ``src`` directory is completely arbitrary and it will be referred in the following as ``msc``.

Within the ``msc`` directory, the user needs to create an ``input`` directory, where the input files and data required to run the MSC model will be contained. The minimum required input file for the MSC model is called ``settings.msm``. Four examples are included in the ``tests`` folder, so that the user can get familiar with running the model. In order to run any of these examples, copy the content of the example directory (e.g. ``rectang_wf_subcritical``) inside the ``input`` directory that has been created within the ``msc`` directory, then run the ``main.m`` in Matlab. 

For additional cases, the user may not be able to generate the ``ABL_averaging.mat`` input file, which contains results from a TOSCA LES of the ABL, nor the ``turbineCurves.m`` file, which contains the turbine-manufacturer curves that allow to run the model in variable-Ct mode. If this is the case, certain parameters can be set inside the ``settings.msm`` file which make these data not required, as they are calculated directly by the MSC model. 

The parameters of the ``settings.msm`` file, together with their description, are summarized in the following table:

| Parameter               | Description                                                                                                 |
|-------------------------|-------------------------------------------------------------------------------------------------------------|
| iterations              | number of coupling iterations (min 1)                                                                       |
| deepArray               | zero: no, two: add deep array effects (only improves in random farm configurations)                         |
| localInduction          | add local turbine induction                                                                                 |
| readTurbineThrust       | thrust from LES (wrapper for TOSCA turbine data)                                                            |
| readTurbineCtPrime      | disk averaged thrust coefficient from LES (wrapper for TOSCA turbine data)                                  |
| periodicSpanwise        | periodize the wind turbine distribution in the spanwise direction                                           |
| ablModel                | one: Nieuwstadt, two: .nc/.dat input from PALM/TOSCA                                                        |
| bkStateName             | only required if ablModel = 2 (PALM input: palm, TOSCA input: tosca)                                        |
| displayMicroScale       | computes also microscale wind around the wind farm (COMPUTATIONALLY EXPENSIVE, for post processing)         |
| wakeModelLogLaw         | zero: use uniform inflow for wake model RHS, one: use log inflow for wake model RHS                         |
| singlePointCoupling     | if set to one performs the original 3LM one-point coupling                                                  |
| mirrorWindFarm          | if set to one mirrors the wind farm w.r.t. the ground                                                       |
| excludeBackground       | if set to one excludes the background wind variation (wake model alone and induction if activated)          |
| rho                     | reference air density [Kg/m3]                                                                               |
| nu                      | kinematic viscosity [m2/s]                                                                                  |
| g                       | mag of gravitational acceleration [m/s2]                                                                    |
| Tref                    | reference ground temperature [K]                                                                            |
| Uref                    | reference wind vector at hub height [m/s]                                                                   |
| Href                    | reference hub height [m]                                                                                    |
| H                       | reference ABL height [m]                                                                                    |
| dTH                     | T jump across capping [K]                                                                                   |
| dTdz                    | lapse rate above capping [K/Km]                                                                             |
| dTdzABL                 | lapse rate below capping [K/Km]                                                                             |
| z0                      | equivalent roughness length [m]                                                                             |
| TI                      | turbulence intensity percent [%]                                                                            |
| lat                     | latitude [deg]                                                                                              |
| xsMeso                  | x mesoscale domain start [Km]                                                                               |
| xeMeso                  | x mesoscale domain end [Km]                                                                                 |
| ysMeso                  | y mesoscale domain start [Km]                                                                               |
| yeMeso                  | y mesoscale domain end [Km]                                                                                 |
| dxMeso                  | mesoscale domain discretization in x dir [m]                                                                |
| dyMeso                  | mesoscale domain discretization in y dir [m]                                                                |
| xsFarm                  | x farm start [Km]                                                                                           |
| ysFarm                  | y farm start [Km]                                                                                           |
| farmInputType           | one: aligned wind farm, two: staggered wind farm, three: read from excel                                    |
| excelName               | name of the excel file where the turbine coords are saved                                                   |
| Sx                      | x turbine spacing in diameters (discarded if farmInputType = three)                                         |
| Sy                      | y turbine spacing in diameters (discarded if farmInputType = three)                                         |
| Nty                     | y number of wind turbines (discarded if farmInputType = three)                                              |
| Ntx                     | x number of wind turbines (discarded if farmInputType = three)                                              |
| D                       | turbine diameter [m] (discarded if farmInputType = three)                                                   |
| zHub                    | turbine hub-height [m] (discarded if farmInputType = three)                                                 |
| useTurbineCurves        | zero: constant Ct and Cp (below), one: read curves from input/turbineCurves.mat (wind, cp, ct tables)       |
| Ct                      | default value of turbine thrust coeff (discarded if farmInputType = three)                                  |
| Cp                      | default value of turbine power coeff (discarded if farmInputType = three)                                   |    

Notably, the ``ABL_averaging.mat`` is required if *ablModel* is set to 2 and *bkStateName* is set to *tosca*. For PALM, different data is required (although it has been a long time that the MSC model has not been tested using the PALM output). If *ablModel* is set to 1, the MSC will run the Nieuwstadt model internally, which provides the required inputs. Notably, the Nieuwstadt model is a very good ABL model that takes both veer and shear into 
account. Using an LES input like PALM or TOSCA does not add any information and it is only suggested when comparing against LES in order to remove the bias produced by different ABL states (and so not strictly due to the MSC model itself). 

Furhtermore, if the user wants to simulate an arbitrary wind farm, an excel file similar to those provided in the ``arbitrary_wf_stable`` and  ``arbitrary_wf_neutral`` examples is required. Please note that the coordinates provided in the file are always rescaled by the MSC model based on the minimum x and y points. Rectancular aligned and staggered wind farms can be defined directly from the ``settings.msm`` file without the need to provide and excel file, upon specifying the number of turbines in x and y, as well as their spacing. 

If the user does not have the turbine curves from the manufacturer (those provided in the examples are for the NREL5MW and they have been derived from LES using TOSCA), a constant Ct can be used by setting ``useTurbineCurves`` to zero. 

## Outputs

The MSC model output is saved into a ``.mat`` file that can be read by Matlab and post-processed with ad-hoc code. The dataset contains the following structs:
- sol   : numerical solution settings and run info
- abl   : atmospheric boundary layer (background state) data 
- meso  : atmopsheric perturbation solution fields 
- farm  : turbine data 
- micro : microscale hub-height fields

Notably, the fields ``micro.u`` and ``micro.v``, calculated by setting ``displayMicroScale`` to 1, are super-expensive to compute (although their computation is parallelized). This is because for each turbine, the wake model and the local induction model (if activated) have to be evaluated for each point of the microscale domain mesh. Note that this is **NOT** how wake models are usually run. In fact, wake models are usually only evaluated at the turbine locations because only the wind speed at these location is required to calculate turbine power. However, the option ``displayMicroScale`` **MAY** be activated if one wishes to create a nice visualization, even if it does not add anything to the MSC model. 

Happy use of the MSC model!
  
