
#knownFluxes

In this folder, you can add the known fluxes. Just make sure that the name of the reaction is consistent with the GEM you will use. 
For example, in this case, we add the oxygen and biomass consumption reactions represented in the iMM904 model as "R_EX_o2_e" and "R_BIOMASS_SC5_notrace".

You must also modify the ```epiflux.py``` in the code line 274 (getAdditionalConstraint(model, v_dic,vo2Known). 
In this section, you must alter or add the "Reaction_ID" of the known fluxes. 
