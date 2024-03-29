import os
import json
import pickle
from icolos.core.workflow_steps.active_learning.surrogate_model import SurrogateModel
from icolos.core.workflow_steps.active_learning.acquisition_function import (
    AcquisitionFunction,
)
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.step_enums import StepBaseEnum, StepGlideEnum
from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.step_enums import (
    StepActiveLearningEnum,
)
from icolos.core.workflow_steps.step import StepBase
from icolos.core.composite_agents.workflow import WorkFlow

from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

_SALE = StepActiveLearningEnum()
_WE = WorkflowEnum()
_SBE = StepBaseEnum
_SGE = StepGlideEnum()


class StepProspectiveREINVENT(StepBase):
    """
    Class to run prospective active learning with REINVENT.
    The scoring function component (e.g., docking via LigPrep + Glide) will be predicted by a surrogate model
    for a fraction of every batch of SMILES generated by REINVENT at a given epoch. The goal is to mitigate
    computational costs of expensive components. The corresponding fraction of SMILES not predicted will be
    sent to the defined oracle and the ground truth labels will be used to augment the surrogate model.
    This process is repeated to constitute an active learning loop
    """

    def __init__(self, **data):
        super().__init__(**data)

    def _run(self, original_smiles: list, save_path: str):
        """
        main method controlling what icolos logic to run out of 3 possibilities:
        1. run normal REINVENT (corresponding to warm-up)
        2. run normal REINVENT with pooling (corresponding to initial pooling epochs)
        3. run active learning REINVENT
        """

        (
            epochs,
            warmup,
            retrain,
            initial_pooling_epochs,
            acquisition_batch_size,
            acquisition_function_name,
            surrogate_model_type,
            oracle_config,
        ) = self._get_run_parameters()

        # extract current epoch
        curr_epoch = self._get_curr_epoch()

        if curr_epoch <= warmup:
            # run "normal" REINVENT, i.e., all generated SMILES are passed to the Oracle
            self._run_normal_reinvent(
                curr_epoch=curr_epoch,
                save_path=save_path,
                original_smiles=original_smiles,
                oracle_config=oracle_config
            )

        elif curr_epoch <= (warmup + initial_pooling_epochs):
            # continue running "normal" REINVENT but start pooling the (SMILES, docking score) pairs
            # to generate the initial training set for the surrogate model
            self._run_normal_reinvent_with_pooling(
                curr_epoch=curr_epoch,
                save_path=save_path,
                original_smiles=original_smiles,
                oracle_config=oracle_config,
            )

        # run active learning REINVENT for the remainder of the experiment
        else:
            self._run_active_learning_reinvent(
                curr_epoch=curr_epoch,
                save_path=save_path,
                original_smiles=original_smiles,
                oracle_config=oracle_config,
                warmup=warmup,
                initial_pooling_epochs=initial_pooling_epochs,
                retrain=retrain,
                surrogate_model_type=surrogate_model_type,
                acquisition_function_name=acquisition_function_name,
                acquisition_batch_size=acquisition_batch_size,
            )

    def _get_morgan_fingerprints(self, smiles) -> list:
        """returns Morgan fingerprints for a list of smiles"""
        molecules = [Chem.MolFromSmiles(s) for s in smiles]
        return [
            Chem.AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048)
            for mol in molecules
        ]

    def _get_curr_epoch(self) -> int:
        """get the current REINVENT epoch number to determine which stage of active learning to execute"""

        global_variables = self.get_workflow_object().header.global_variables
        curr_epoch = int(global_variables["step_id"])
        # curr_epoch + 1 to start at epoch 1
        return curr_epoch + 1

    def _save_state(
        self,
        curr_epoch: int,
        save_path: str,
        original_smiles: list,
        fingerprints: list,
        labels: list,
        pooled_smiles=None,
        pooled_labels=None,
        pooled_fingerprints=None,
        surrogate=None,
        save_pool=False,
        non_acquired_smiles=None,
        non_acquired_labels=None,
    ):
        """stores intermediate results after every epoch"""
        # save the pooled compounds to be used as training data for the surrogate model
        if save_pool:
            # make sure pooled_compounds are smiles, but we should also save the enumerated forms too?
            pooled_smiles.extend(original_smiles)
            pooled_labels.extend(labels)
            df = pd.DataFrame({"smiles": pooled_smiles, "labels": pooled_labels})
            df.to_csv(os.path.join(save_path, f"pooled_data/epoch_{curr_epoch}.csv"))

            # if using PI or EI, save the current best
            af = self._get_additional_setting(_SALE.ACQUISITION_FUNCTION, default=_SALE.RANDOM)
            if (af == _SALE.PI) or (af == _SALE.EI):
                current_best = min(pooled_labels)
                with open("current_best.txt", "w+") as f:
                    f.write(f"{str(current_best)}\n{str(curr_epoch)}")

            pooled_fingerprints.extend(fingerprints)
            pkl_path = os.path.join(save_path, f"pooled_data/epoch_{curr_epoch}_fingerprints.sav")
            pickle.dump(pooled_fingerprints, open(pkl_path, "wb"))

        # save the non-acquired SMILES and the predicted labels, if applicable
        if non_acquired_smiles is not None:
            df = pd.DataFrame(
                {"smiles": non_acquired_smiles, "labels": non_acquired_labels}
            )
            df.to_csv(
                os.path.join(save_path, f"predicted_smiles/epoch_{curr_epoch}.csv")
            )

        # save the surrogate model state
        if isinstance(surrogate, SurrogateModel):
            pkl_path = os.path.join(
                save_path, f"surrogate_states/epoch_{curr_epoch}.sav"
            )
            pickle.dump(surrogate, open(pkl_path, "wb"))

    def _read_state(self, curr_epoch, save_path: str, train_negative_pool: False):
        """extract the previous epoch's intermediate results to be used in the current epoch"""
        # read the saved state CSV file and extract the pooled SMILES and labels
        # at the first pooling epoch, there will be neither pooled compounds nor surrogate state
        try:
            # curr_epoch - 1 to ensure the previous state prior to the current active learning iteration is read
            df = pd.read_csv(
                os.path.join(save_path, f"pooled_data/epoch_{curr_epoch-1}.csv")
            )
            pooled_smiles, pooled_labels = list(df["smiles"]), list(df["labels"])

            # ------------------------------------------------------------------------------------
            # TODO: temporary implementation to train on warm-up data
            if train_negative_pool:
                files = os.listdir("docking_scores")
                files = sorted(files, key=lambda x: int(x.split("_")[0]))
                for file in files:
                    df = pd.read_csv(os.path.join("docking_scores", file))
                    pooled_smiles.extend(df["original_smiles"])
                    pooled_labels.extend(df["docking_score"])
                self._logger.log("added oracle scores", _LE.DEBUG)
                pooled_fingerprints = self._get_morgan_fingerprints(smiles=pooled_smiles)
            # --------------------------------------------------------------------------------------
            else:
                pkl_path = os.path.join(save_path, f"pooled_data/epoch_{curr_epoch-1}_fingerprints.sav")
                pooled_fingerprints = pickle.load(open(pkl_path, 'rb'))
                # TODO: could replace these magic strings

        except Exception:
            pooled_smiles, pooled_fingerprints, pooled_labels = [], [], []

        try:
            surrogate = pickle.load(
                open(
                    os.path.join(
                        save_path, f"surrogate_states/epoch_{curr_epoch-1}.sav"
                    ),
                    "rb",
                )
            )
        except Exception:
            surrogate = None

        return pooled_smiles, pooled_fingerprints, pooled_labels, surrogate

    def _query_oracle(self, save_path: str, original_smiles, oracle_config) -> list:
        """executes the oracle Workflow object and returns the scores"""

        oracle_workflow = self._construct_oracle_workflow(oracle_config=oracle_config)
        # copy and save the REINVENT's batch of SMILES to a temporary save folder
        self._create_manipulable_compounds_file(
            original_smiles=original_smiles, save_path=save_path
        )
        oracle_workflow.execute()

        # extract the best scores per compound and its enumerations
        scores = self._extract_scores_from_oracle(
            original_smiles=original_smiles, oracle_workflow=oracle_workflow
        )

        return scores

    def _create_manipulable_compounds_file(self, original_smiles, save_path: str):
        """
        this method copies the compounds generated by REINVENT into a new directory.
        The compounds here can be freely manipulated, in particular, partitioning into
        compounds to acquire and compounds to predict via the surrogate model
        """

        path = os.path.join(save_path, "reinvent_batches")
        with open(os.path.join(path, "oracle_compounds.smi"), "w+") as f:
            for smiles in original_smiles:
                f.write(f"{smiles}\n")

    def _extract_scores_from_oracle(
        self, original_smiles, oracle_workflow: WorkFlow
    ) -> list:
        """
        extracts the "best" score from the oracle and is particularly relevant in cases where more than 1 score
        is given to a SMILES. This is usually taken care of by icolos write-out and is done here only because
        the training data for the surrogate model (SMILES, score pairs) needs to be stored
        """
        scores = []
        for smiles, compound in zip(
            original_smiles, oracle_workflow.get_steps()[-1].get_compounds()
        ):
            # TODO: assumes lower the score the better. Add a parameter in the JSON to control this for generalizability
            try:
                    scores_tracker = []
                    for enumeration in compound.get_enumerations():
                        for conformer in enumeration.get_conformers():
                            # docking scores here are with Epik corrections
                            # TODO: docking score is hard-coded at the moment, need an icolos parameter to specifiy this
                            score = float(
                                conformer.get_molecule().GetProp(_SGE.GLIDE_DOCKING_SCORE)
                            )
                            scores_tracker.append(score)
                    # TODO: this may cause an error if list is empty --> enumeration not added in dummy conf?
                    #  try block handles this for the time being
                    scores.append(min(scores_tracker))
            except Exception:
                scores.append(float(0.0))

        return scores

    def _write_reinvent_feedback(self, original_smiles, labels):
        """this method manually writes the REINVENT feedback"""
        path = self.get_workflow_object().header.global_variables["output_json_path"]

        # below is the form expected by REINVENT
        # extract the oracle label
        values_key = self._get_additional_setting(_SALE.ORACLE_LABEL)
        try:
            output_dict = {
                "results": [{"values_key": values_key, "values": labels}],
                "names": [idx for idx in range(len(original_smiles))],
            }
        # if for some reason, the number of SMILES is not equal to the number of labels, labels (0.0 float)
        # are manually appended to prevent code error. This is a temporary fix for the initial experiments
        # and one should ensure in the icolos step that each SMILES is accounted for
        except Exception:
            missing_values = len(original_smiles) - len(values_key)
            append_values = ["dummy SMILES" for idx in range(missing_values)]
            values_key.extend(append_values)

            missing_labels = len(original_smiles) - len(labels)
            append_labels = [0.0 for idx in range(missing_labels)]
            labels.extend(append_labels)

            output_dict = {
                "results": [{"values_key": values_key, "values": labels}],
                "names": [idx for idx in range(len(original_smiles))],
            }

            self._logger.log("REINVENT Feedback Missing Values. Manually Added.", _LE.DEBUG)

        with open(path, "w+") as f:
            json.dump(output_dict, f, indent=2)

    def _run_normal_reinvent(
        self, curr_epoch, save_path, original_smiles, oracle_config
    ):
        """runs a normal REINVENT epoch"""
        scores = self._query_oracle(
            save_path=save_path,
            original_smiles=original_smiles,
            oracle_config=oracle_config,
        )

        self._save_state(
            curr_epoch=curr_epoch,
            original_smiles=original_smiles,
            fingerprints=None,
            labels=scores,
            save_path=save_path,
        )
        self._write_reinvent_feedback(original_smiles=original_smiles, labels=scores)

    def _run_normal_reinvent_with_pooling(
        self, curr_epoch, save_path: str, original_smiles, oracle_config
    ):
        """runs a normal REINVENT epoch and stores all (SMILES, scores) pairs as training data for the surrogate"""
        # read the current state
        pooled_smiles, pooled_fingerprints, pooled_labels, surrogate = self._read_state(
            curr_epoch=curr_epoch, save_path=save_path, train_negative_pool=False
        )

        scores = self._query_oracle(
            save_path=save_path,
            original_smiles=original_smiles,
            oracle_config=oracle_config,
        )

        # generate the Morgan fingerprints for saving
        fingerprints = self._get_morgan_fingerprints(smiles=original_smiles)

        self._save_state(
            curr_epoch=curr_epoch,
            original_smiles=original_smiles,
            fingerprints=fingerprints,
            labels=scores,
            pooled_smiles=pooled_smiles,
            pooled_labels=pooled_labels,
            pooled_fingerprints=pooled_fingerprints,
            save_pool=True,
            save_path=save_path,
        )

        self._write_reinvent_feedback(original_smiles=original_smiles, labels=scores)

    def _run_active_learning_reinvent(
        self,
        curr_epoch,
        save_path: str,
        original_smiles,
        oracle_config,
        warmup,
        initial_pooling_epochs,
        retrain,
        surrogate_model_type,
        acquisition_function_name,
        acquisition_batch_size,
    ):
        """runs active learning REINVENT"""

        # the if-else block below was a quick implementation of training the initial surrogate model state
        # using the "negative" data during warm-up. The train_negative_pool boolean of _read_state controls this.
        # In the current code state, it is set to False in both the if-else blocks but could be set to True in
        # the if block to train with negative data
        # TODO: should have a boolean parameter in the icolos json specify whether or not to train with negative data
        if curr_epoch == (warmup + initial_pooling_epochs + 1):
            pooled_smiles, pooled_fingerprints, pooled_labels, surrogate = self._read_state(
                curr_epoch=curr_epoch, save_path=save_path, train_negative_pool=False
            )
        else:
            pooled_smiles, pooled_fingerprints, pooled_labels, surrogate = self._read_state(
                curr_epoch=curr_epoch, save_path=save_path, train_negative_pool=False
            )

        # model training is controlled by the "retrain" parameter. The only exception to this is when
        # curr_epoch = (warmup + initial_pooling_epochs + 1), i.e., the epoch right after initial pooling
        # this epoch corresponds to the initial training of the surrogate model and is also the first epoch
        # where we begin saving the pickled surrogate model
        if (
            curr_epoch == (warmup + initial_pooling_epochs + 1)
            or curr_epoch % retrain == 0
        ):
            surrogate = SurrogateModel(model_type=surrogate_model_type)
            surrogate.fit(pooled_fingerprints, pooled_labels)
        # at this point, whether the surrogate was retrained or not, partition the compounds into acquired and
        # non-acquired. Send acquired to oracle, predict the rest, and send the concatenated data back to REINVENT
        acquisition_function = AcquisitionFunction(
            surrogate=surrogate,
            function=acquisition_function_name,
            acquisition_batch_size=acquisition_batch_size,
        )

        # select points to acquire based on acquisition function
        acquired_smiles, non_acquired_smiles = acquisition_function.partition_compounds(
            original_smiles=original_smiles, fingerprints_pool=pooled_fingerprints
        )

        # send the acquired compounds to the oracle
        acquired_labels = self._query_oracle(
            save_path=save_path,
            original_smiles=acquired_smiles,
            oracle_config=oracle_config,
        )
        # predict the non-acquired compounds and cast to type list from np.array for compatibility with other methods
        non_acquired_labels = list(
            surrogate.predict(self._get_morgan_fingerprints(smiles=non_acquired_smiles))
        )

        # generate the Morgan fingerprints of the acquired compounds for saving
        acquired_fingerprints = self._get_morgan_fingerprints(smiles=acquired_smiles)

        # save the current state for the next iteration
        # note: original_smiles and labels here correspond to the ***acquired*** data.
        # Concatenating the acquired data to the current pool is handled by save_state.
        # The non-acquired data is also saved for analysis at the end of the experiment
        self._save_state(
            curr_epoch=curr_epoch,
            save_path=save_path,
            original_smiles=acquired_smiles,
            fingerprints=acquired_fingerprints,
            labels=acquired_labels,
            pooled_smiles=pooled_smiles,
            pooled_labels=pooled_labels,
            pooled_fingerprints=pooled_fingerprints,
            non_acquired_smiles=non_acquired_smiles,
            non_acquired_labels=non_acquired_labels,
            surrogate=surrogate,
            save_pool=True,
        )

        # concatenate the acquired and non-acquired points to provide REINVENT feedback
        # from REINVENT's perspective, it has no notion of acquired vs. predicted. It simply receives scores back
        concatenated_smiles = acquired_smiles + non_acquired_smiles
        concatenated_labels = acquired_labels + non_acquired_labels
        self._write_reinvent_feedback(
            original_smiles=concatenated_smiles, labels=concatenated_labels
        )

    def _get_run_parameters(self):
        """returns all the experiment parameters"""
        # default 1000 REINVENT epochs if running standard reinforcement learning for small molecules
        epochs = self._get_additional_setting(_SALE.EPOCHS, default=1000)
        warmup = self._get_additional_setting(_SALE.WARMUP, default=200)
        retrain = self._get_additional_setting(_SALE.RETRAIN, default=5)
        initial_pooling_epochs = self._get_additional_setting(
            _SALE.INITIAL_POOLING_EPOCHS, default=5
        )
        acquisition_batch_size = self._get_additional_setting(
            _SALE.ACQUISITION_BATCH_SIZE, default=40
        )
        acquisition_function = self._get_additional_setting(
            _SALE.ACQUISITION_FUNCTION, default=_SALE.RANDOM
        )
        surrogate_model_type = self._get_additional_setting(
            _SALE.SURROGATE_MODEL_TYPE, default=_SALE.RANDOM_FOREST_REGRESSOR
        )
        oracle_config = self._get_additional_setting(_SALE.ORACLE_CONFIG)

        return (
            epochs,
            warmup,
            retrain,
            initial_pooling_epochs,
            acquisition_batch_size,
            acquisition_function,
            surrogate_model_type,
            oracle_config,
        )

    def _construct_oracle_workflow(self, oracle_config) -> WorkFlow:
        """constructs an icolos Workflow object for the oracle (essentially a sub-Workflow)"""

        # load oracle configuration
        with open(oracle_config, "r") as file:
            conf = file.read().replace("\r", "").replace("\n", "")
            conf = json.loads(conf)

        # generate oracle workflow object
        oracle_workflow = WorkFlow(**conf[_WE.WORKFLOW])
        # copy over the header from the parent workflow to extract the icolos global variables
        header = self.get_workflow_object().header
        oracle_workflow.header = header
        oracle_workflow.initialize()

        return oracle_workflow

    def _create_save_folders(self, save_path: str):
        """creates necessary directories for storage of intermediate results"""
        if not os.path.exists(save_path):
            os.mkdir(save_path)

        if not os.path.exists(os.path.join(save_path, "pooled_data")):
            os.mkdir(os.path.join(save_path, "pooled_data"))

        if not os.path.exists(os.path.join(save_path, "predicted_smiles")):
            os.mkdir(os.path.join(save_path, "predicted_smiles"))

        if not os.path.exists(os.path.join(save_path, "surrogate_states")):
            os.mkdir(os.path.join(save_path, "surrogate_states"))

        if not os.path.exists(os.path.join(save_path, "reinvent_batches")):
            os.mkdir(os.path.join(save_path, "reinvent_batches"))

    def execute(self):
        # extract the original SMILES from the REINVENT batch
        original_smiles = [
            compound.get_enumerations()[0].get_original_smile()
            for compound in self.data.compounds
        ]
        # extract the save directory specified to store intermediate results
        save_path = self._get_additional_setting(_SALE.SAVE_DIR)
        # create save folders to store pooled SMILES, pooled fingerprints, and surrogate model states
        self._create_save_folders(save_path=save_path)

        self._run(original_smiles=original_smiles, save_path=save_path)
