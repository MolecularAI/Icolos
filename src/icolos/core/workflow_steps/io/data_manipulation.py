from typing import List, Union
from pydantic import BaseModel

from icolos.core.containers.compound import unroll_conformers
from icolos.core.step_utils.structcat_util import StructcatUtil
from icolos.core.step_utils.structconvert import StructConvert
from icolos.utils.enums.program_parameters import (
    OpenBabelEnum,
    SchrodingerExecutablesEnum,
)
from icolos.utils.enums.step_enums import (
    StepDataManipulationEnum,
    StepBaseEnum,
    StepFilterEnum,
)
from icolos.core.workflow_steps.io.base import StepIOBase
import os
from icolos.core.workflow_steps.step import _LE
import numpy as np

_SBE = StepBaseEnum
_SDM = StepDataManipulationEnum()
_SEE = SchrodingerExecutablesEnum()
_OE = OpenBabelEnum()
_SFE = StepFilterEnum()


class StepDataManipulation(StepIOBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        # extend parameters
        if _SDM.ACTION not in self.settings.additional.keys():
            self.settings.additional[
                _SDM.ACTION
            ] = _SDM.ACTION_ATTACH_CONFORMERS_AS_EXTRA
            self._logger.log(
                f"Action not specified, defaulting to {_SDM.ACTION_ATTACH_CONFORMERS_AS_EXTRA}.",
                _LE.WARNING,
            )

    def _attach_conformers_as_extra(self):
        # load data to match from previous step (note: no other input supported here, to avoid redudancy
        # with standard input preparation)
        match_compounds = (
            self.get_workflow_object()
            .find_step_by_step_id(self.settings.additional[_SDM.MATCH_SOURCE])
            .clone_compounds()
        )

        # unroll for convenience, attach matches to input conformers as extra data
        match_conformers = unroll_conformers(match_compounds)
        for comp in self.get_compounds():
            for enum in comp:
                for conf in enum:
                    list_matched = [
                        c
                        for c in match_conformers
                        if conf.get_index_string() == c.get_index_string()
                    ]
                    conf.add_extra_data(key=_SDM.KEY_MATCHED, data=list_matched)
                    self._logger.log(
                        f"Added {len(list_matched)} conformers as extra data to conformer {conf.get_index_string()}.",
                        _LE.DEBUG,
                    )

    def _convert_mae_to_pdb(self):
        converter = StructConvert(prefix_execution=_SEE.SCHRODINGER_MODULE)
        tmp_dir = self._make_tmpdir()

        # find the mae files from the input step and convert to pdb
        for file in self.data.generic.get_files_by_extension("mae"):
            file.write(tmp_dir)
            output_file = file.get_file_name().split(".")[0] + ".pdb"
            converter.mae2pdb(
                os.path.join(tmp_dir, file.get_file_name()),
                os.path.join(tmp_dir, output_file),
            )
        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)

    def _assemble_complexes(self):
        concatenator = StructcatUtil(
            prefix_execution=_SEE.SCHRODINGER_MODULE, backend=_OE.OBABEL
        )
        assert os.path.isfile(self.settings.additional[_SDM.RECEPTOR])
        # create a tmpdir to work in
        tmp_dir = self._make_tmpdir()
        # get compounds from previous step
        conformers = self._unroll_compounds(self.get_compounds(), level="conformers")
        for conf in conformers:
            path = os.path.join(tmp_dir, f"{conf.get_index_string()}.sdf")
            mol = conf.get_molecule()
            conf.write(path)
            concatenator.concatenate(
                input_files=[
                    self.settings.additional[_SDM.RECEPTOR],
                    path,
                ],
                output_file=os.path.join(tmp_dir, f"{conf.get_index_string()}.pdb"),
            )
        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)

    def _filter_compounds(self):
        """modifies set of input compounds according to the specification provided in the config block"""
        # TODO: support ranking structures based on generic data

        top_n = self.settings.additional[_SFE.RETURN_N]
        reverse = self.settings.additional[_SFE.HIGHEST_IS_BEST]
        criteria = (
            self.settings.additional[_SFE.CRITERIA]
            if _SFE.CRITERIA in self.settings.additional.keys()
            else None
        )
        aggregation = (
            self.settings.additional[_SFE.AGGREGATION]
            if _SFE.AGGREGATION in self.settings.additional.keys()
            else "sum"
        )

        top_conformer_list = []
        for compound in self.data.compounds:
            # filter by enumeration first - return a list of the top scoring conformers for that enumeration
            # this is the normal running mode, as opposed to sorting by compound, regardless of the enumeration it came from
            for enumeration in compound.get_enumerations():
                enumeration.sort_conformers(
                    by_tag=criteria, reverse=reverse, aggregation=aggregation
                )
                top_confs = enumeration.get_conformers()[:top_n]
                enumeration.clear_conformers()
                enumeration.add_conformers(top_confs)
                # replace that enumeration's conformers with the sorted list.
                # if filtering at conformer level i.e. regardless of enumeration
                if self.settings.additional[_SFE.FILTER_LEVEL] == _SFE.COMPOUNDS:
                    for conf in top_confs:
                        top_conformer_list.append(conf)
        if self.settings.additional[_SFE.FILTER_LEVEL] == _SFE.COMPOUNDS:
            # sort the top conformers from each enumeration and attach the top n conformers to their respective enumeration, get rid of the rest
            # sorted_top_confs = sorted(top_conformer_list,
            #                           key=lambda x: x.get_molecule().GetProp(self.settings.additional[_SFE.CRITERIA]),
            #                           reverse=reverse)[:top_n]
            # sort conformers
            sorted_top_confs = self._sort_conformers(
                conformers=top_conformer_list,
                by_tag=criteria,
                reverse=reverse,
                aggregation=aggregation,
            )
            for compound in self.data.compounds:
                for enum in compound.get_enumerations():
                    enum.clear_conformers()
            for conf in sorted_top_confs:
                enum = conf.get_enumeration_object()
                enum.add_conformer(conf)

    def _sort_conformers(
        self,
        conformers,
        by_tag: Union[str, List[str]],
        reverse: bool = True,
        aggregation="sum",
    ):
        if isinstance(by_tag, list) and len(by_tag) == 1:
            by_tag = by_tag[0]

        if isinstance(by_tag, str):
            # sorting according to a single tag
            conformers = sorted(
                conformers,
                key=lambda x: float(x.get_molecule().GetProp(by_tag)),
                reverse=reverse,
            )
            return conformers
            # self._conformers = conformers
            # self.reset_conformer_ids()
        elif isinstance(by_tag, list):
            # need to normalise the values, calculate max and min of each tag for that series of conformers provided
            # this would allow us to compare across a series, i.e. scoring and ranking the output of all conformers in an enumeration from Glide
            def normalise_tag(value, tag):
                all_tag_values = [
                    float(conf.get_molecule().GetProp(tag)) for conf in conformers
                ]
                if len(all_tag_values) == 1:
                    return value
                else:

                    max_tag = np.max(all_tag_values)
                    min_tag = np.min(all_tag_values)
                    return (float(value) - min_tag) / (max_tag - min_tag)

            # if we specify multiple tags, aggregate according the the provided aggregation function
            if aggregation == "sum":
                # sort by the sum of the normalised tags,
                conformers = sorted(
                    conformers,
                    key=lambda x: np.sum(
                        [
                            float(normalise_tag(x.get_molecule().GetProp(i), i))
                            for i in by_tag
                        ]
                    ),
                    reverse=reverse,
                )
                return conformers
            elif aggregation == "product":
                conformers = sorted(
                    conformers,
                    key=lambda x: np.product(
                        [
                            float(normalise_tag(x.get_molecule().GetProp(i), i))
                            for i in by_tag
                        ]
                    ),
                    reverse=reverse,
                )
                return conformers
            else:
                raise AttributeError(
                    "Only sum or product aggregation modes are currently supported - ABORT"
                )

    def execute(self):
        if (
            self.settings.additional[_SDM.ACTION]
            == _SDM.ACTION_ATTACH_CONFORMERS_AS_EXTRA
        ):
            self._attach_conformers_as_extra()
        elif self.settings.additional[_SDM.ACTION] == _SDM.ACTION_NO_ACTION:
            n_comp, n_enum, n_conf = self.get_compound_stats()
            self._logger.log(
                f'Data manipulation step type "no_action" for {n_comp} compounds with {n_enum} enumerations with {n_conf} conformers completed.',
                _LE.INFO,
            )
        elif self.settings.additional[_SDM.ACTION] == _SDM.CONVERT_MAE_TO_PDB:
            self._convert_mae_to_pdb()
        elif self.settings.additional[_SDM.ACTION] == _SDM.ASSEMBLE_COMPLEXES:
            # take pose conformers (sd format) and concatenate with pdb file
            self._assemble_complexes()
        elif self.settings.additional[_SDM.ACTION] == _SDM.COLLECT_ITERATOR_RESULTS:
            # average the results coming from all iterations of the step
            raise NotImplementedError
        elif self.settings.additional[_SDM.ACTION] == _SDM.FILTER:
            self._filter_compounds()
        else:
            raise ValueError(
                f'Action "{self.settings.additional[_SDM.ACTION]}" not supported.'
            )
