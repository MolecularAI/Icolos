from copy import deepcopy

import numpy as np
from typing import List

from pydantic import BaseModel

from icolos.core.containers.compound import Enumeration, Conformer

from icolos.utils.enums.step_enums import StepBoltzmannWeightingEnum
from icolos.core.workflow_steps.calculation.base import StepCalculationBase

from icolos.utils.general.convenience_functions import *
from icolos.utils.constants import *

_SBWE = StepBoltzmannWeightingEnum()


class StepBoltzmannWeighting(StepCalculationBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

    def _get_relative_energy_values(
        self, conformers: List[Conformer], property_name: str
    ) -> List[float]:
        values = [float(c.get_molecule().GetProp(property_name)) for c in conformers]
        min_val = min(values)
        relative_values = [value - min_val for value in values]
        return relative_values

    def _obtain_factors(self, relative_values: List[float]) -> List[float]:
        # calculate individual Boltzmann factors
        individual_factors = [
            np.exp((-1 * val / (CONSTANT_KB * CONSTANT_T))) for val in relative_values
        ]

        # calculate and return Boltzmann factors
        sum_factors = sum(individual_factors)
        factors = [val / sum_factors for val in individual_factors]
        return factors

    def _calculate_Boltzmann_factors(
        self, enumeration: Enumeration, parameters: dict
    ) -> List[str]:
        list_properties = parameters[_SBWE.PROPERTIES]
        list_output_names = []
        for prop in list_properties:
            # (1) get the relative values for this property (e.g. solvent) for all conformers in respect to the one
            # with the minimal energy
            relative_prop_values = self._get_relative_energy_values(
                conformers=enumeration.get_conformers(),
                property_name=prop[_SBWE.PROPERTIES_INPUT],
            )

            # (2) calculate the Boltzmann factors for this property
            boltzmann_factors = self._obtain_factors(
                relative_values=relative_prop_values
            )

            # (3) add the Boltzmann factors to the conformers as a tag
            for c, bm_factor in zip(enumeration.get_conformers(), boltzmann_factors):
                c.get_molecule().SetProp(prop[_SBWE.PROPERTIES_OUTPUT], str(bm_factor))
            list_output_names.append(prop[_SBWE.PROPERTIES_OUTPUT])
        return list_output_names

    def _do_Boltzmann_weighting(self, conformers: List[Conformer], weightings: dict):
        input_tags = weightings[_SBWE.WEIGHT_INPUT]
        output_prefix = nested_get(
            weightings, _SBWE.WEIGHT_OUTPUT_PREFIX, default="bf_weighted"
        )
        properties = weightings[_SBWE.WEIGHT_PROPERTIES]
        for prop in properties:
            for inp_tag in input_tags:
                new_tag_name = "_".join([output_prefix, inp_tag, prop])
                products = []
                for conformer in conformers:
                    conf = conformer.get_molecule()
                    products.append(
                        float(conf.GetProp(prop)) * float(conf.GetProp(inp_tag))
                    )
                for conformer in conformers:
                    conformer.get_molecule().SetProp(new_tag_name, str(sum(products)))

    def execute(self):
        parameters = deepcopy(self.settings.arguments.parameters)
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                if self._input_object_empty(enumeration):
                    continue

                # get the name of the Boltzmann properties / solvents and annotate the factors
                _ = self._calculate_Boltzmann_factors(enumeration, parameters)

                # for each property and each weighting, add the respective tags
                self._do_Boltzmann_weighting(
                    conformers=enumeration.get_conformers(),
                    weightings=parameters[_SBWE.WEIGHT],
                )
