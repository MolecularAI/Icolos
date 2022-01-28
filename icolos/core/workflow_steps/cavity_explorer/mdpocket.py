from icolos.utils.general.parallelization import Parallelizer, SubtaskContainer
from pydantic import BaseModel
from icolos.core.workflow_steps.cavity_explorer.base import StepCavityExplorerBase
from icolos.utils.enums.step_enums import StepCavExploreEnum
from icolos.utils.execute_external.execute import Executor
from icolos.core.workflow_steps.step import _LE
from sklearn.cluster import DBSCAN
from collections import Counter
import numpy as np
import re
import os

_SFP = StepCavExploreEnum()


class StepMDpocket(StepCavityExplorerBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        # self._initialize_backend(executor=MPIExecutor)
        self._initialize_backend(executor=Executor)

        # set max_length_sublists to 1
        self.execution.parallelization.max_length_sublists = 1

    def _create_density_grid_file(self, tmp_dir: str, iso_value: float):
        """creates a density grid from the .dx-file into a .pdb-file, heavily influenced by extractISOPdb.py provided
        by fpocket"""
        density_file = [
            file for file in os.listdir(tmp_dir) if file.endswith("dens_grid.dx")
        ]
        assert len(density_file) == 1
        density_file = density_file[0]

        outfile = os.path.join(tmp_dir, f"iso{iso_value}.pdb")

        with open(os.path.join(tmp_dir, density_file), "r") as f:
            # get the axis that shows the most variation during the trajectory, this will be the leading axis
            # read the header - here is an example
            header = ""
            tmp = f.readline()
            while tmp[0] != "o":
                header = header + tmp
                tmp = f.readline()

            # read the grid size
            r = re.compile("\w+")
            gsize = r.findall(tmp)
            gsize = [int(gsize[-3]), int(gsize[-2]), int(gsize[-1])]

            # read the origin of the system
            line = f.readline().split()
            origin = [float(line[-3]), float(line[-2]), float(line[-1])]

            # read grid space
            line = f.readline().split()
            deltax = [float(line[-3]), float(line[-2]), float(line[-1])]
            line = f.readline().split()
            deltay = [float(line[-3]), float(line[-2]), float(line[-1])]
            line = f.readline().split()
            deltaz = [float(line[-3]), float(line[-2]), float(line[-1])]

            # pay attention here, this assumes always orthogonal normalized space, but normally it should be ok
            delta = np.array([deltax[0], deltay[1], deltaz[2]])

            # read the number of data
            f.readline()
            r = re.compile("\d+")
            n_entries = int(r.findall(f.readline())[2])

            if n_entries != gsize[0] * gsize[1] * gsize[2]:
                raise AssertionError(
                    "Error reading the file. The number of expected data points does not correspond to the number of "
                    "labeled data points in the header."
                )
            # create a 3D numpy array filled up with 0
            # initiate xyz counter for reading the grid data
            z = 0
            y = 0
            x = 0

            self._logger.log("Reading grid file...", _LE.DEBUG)

            with open(outfile, "w") as f_out:
                counter = 1
                for _ in range(n_entries // 3):
                    c = f.readline().split()
                    if len(c) != 3:
                        self._logger.log("error reading grid data", _LE.ERROR)
                        raise AssertionError
                    for i in range(3):
                        if (0 > iso_value > float(c[i])) or (
                            0 < iso_value < float(c[i])
                        ):
                            # f_out.write(f"ATOM  {counter}  C   PTH     1   {origin[0] + float(x) * delta[0]} {origin[1] + float(y) * delta[1]} {origin[2] + float(z) * delta[2]} 0.00 0.00\n")
                            f_out.write(
                                "ATOM  %5d  C   PTH     1    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                                % (
                                    counter,
                                    origin[0] + float(x) * delta[0],
                                    origin[1] + float(y) * delta[1],
                                    origin[2] + float(z) * delta[2],
                                    0.0,
                                    0.0,
                                )
                            )
                            counter += 1
                        z += 1
                        if z >= gsize[2]:
                            z = 0
                            y += 1
                            if y >= gsize[1]:
                                y = 0
                                x += 1

        self._logger.log(f"Finished writing {outfile}", _LE.DEBUG)

    def _cluster_pockets(self, tmp_dir, eps, min_samples, threshold, iso_value):
        """
        Clusters points from the initial MDpocket density grid, at a certain iso value
        """
        iso_file = os.path.join(tmp_dir, f"iso{iso_value}.pdb")
        with open(iso_file, "r") as f:
            # collects the data from the pdb-file (x,y,z coordinates)
            data = {
                (line[5:11].strip()): (
                    line[30:38].strip(),
                    line[38:46].strip(),
                    line[46:54].strip(),
                )
                for line in f.readlines()
            }
        db = DBSCAN(eps=eps, min_samples=min_samples).fit(
            np.array(list(data.values())).astype(np.float64)
        )

        labels = db.labels_
        data_ = np.array(list(data.values())).astype(np.float64)
        db.fit_predict(data_)

        self._logger.log(
            f"Number of clusters found for eps = {eps}, iso = {iso_value}, min_samples = {min_samples} and threshold = {threshold} is: {len(set(db.labels_))}",
            _LE.DEBUG,
        )

        pockets_report = Counter(db.labels_)
        filtered_pockets = []
        filtered_data = {}
        filtered_labels = []

        # sorts out the pockets with more than threshold points
        for k, v in pockets_report.items():
            if v > self.threshold and k >= 0:
                filtered_pockets.append(k)

        # get the keys and labels for each data point
        res = list(zip(list(data.keys()), labels))

        # get lists with the data and labels for the filtered pockets
        for pocket in filtered_pockets:
            for (index, label) in res:
                if label == pocket:
                    filtered_data[index] = data.get(index)
                    filtered_labels.append(label)

        self._logger.log(
            f"PocketIDs having more than {self.threshold} points are: {filtered_pockets}",
            _LE.DEBUG,
        )
        self._logger.log(
            f"The number of filtered pockets is: {len(filtered_pockets)}", _LE.DEBUG
        )
        return data, labels, filtered_data, filtered_labels, pockets_report

    def _save_pocket_files(self, tmp_dir, data, labels):
        """saves the individual pockets as individual pdbs to be used with mdpocket"""
        iso_file = os.path.join(tmp_dir, f"iso{self.iso_value}.pdb")

        # define labels and indices
        res = list(zip(list(data.keys()), labels))
        with open(iso_file, "r") as f:
            original_lines = f.readlines()
        # filter out the indices of all pockets except outliers
        indices = list(set([l for l in labels if l >= 0]))

        # save the pocket-pdbs - these are passed with --selected-pocket arg later
        for label in indices:
            with open(os.path.join(tmp_dir, f"pocket_{label}.pdb"), "w") as f:
                for (index, lab) in res:
                    if lab == label:
                        f.write(original_lines[int(index) - 1])

    def _run_mdpocket_selected_pocket(self, tmp_dir):
        """runs the second mdpocket command for fpocket3"""
        pocket_files = [
            file
            for file in os.listdir(tmp_dir)
            if file.endswith(".pdb") and "pocket_" in file
        ]
        argument_dicts = []
        for file in pocket_files:
            arguments = self._parse_arguments(
                flag_dict={
                    "--trajectory_file": os.path.join(
                        tmp_dir,
                        self.data.generic.get_argument_by_extension(self.format_),
                    ),
                    "--trajectory_format": self.format_,
                    "--selected_pocket": os.path.join(
                        tmp_dir, os.path.join(tmp_dir, file)
                    ),
                    "-f": self.data.generic.get_argument_by_extension("pdb"),
                    "-o": file.split(".")[0],
                }
            )
            argument_dicts.append(arguments)

        fpocket_parallelizer = Parallelizer(func=self._execute_mdpocket)
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(argument_dicts)

        while self._subtask_container.done() is False:
            next_batch = self._get_sublists(get_first_n_lists=self._get_number_cores())

            batch_dirs, batch_args = self._prepare_batch_inputs(next_batch, tmp_dir)

            fpocket_parallelizer.execute_parallel(
                tmp_dir=batch_dirs, arguments=batch_args
            )
            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]

    def _prepare_batch_inputs(self, batch, tmp_dir):
        tmp_dirs = []
        args = []
        for next_subtask_list in batch:
            tmp_dirs.append(tmp_dir)
            for (
                subtask
            ) in (
                next_subtask_list
            ):  # enforced only one task per subtask, otherwise it makes no sense
                args.append(subtask.data)  # append the arguments list
        return tmp_dirs, args

    def _execute_mdpocket(self, tmp_dir, arguments):

        self._backend_executor.execute(
            command=_SFP.MDPOCKET_COMMAND,
            arguments=arguments,
            location=tmp_dir,
            check=True,
        )

    def execute(self):

        tmp_dir = self._make_tmpdir()
        # print(paths)
        self._write_input_files(tmp_dir)
        # set some constants from the arguments
        self._set_mdpocket_args()

        # execute the initial mdpocket job (without a specific pocket) to produce the .dx file
        mdpocket_run1_args = self._parse_arguments(
            flag_dict={
                "--trajectory_file": os.path.join(
                    tmp_dir, self.data.generic.get_argument_by_extension(self.format_)
                ),
                "--trajectory_format": self.format_,
                "-f": os.path.join(
                    tmp_dir, self.data.generic.get_argument_by_extension("pdb")
                ),
            }
        )

        # run the first command, produce the dx file and a bunch of pocket_n.pdb pocket topology files
        self._execute_mdpocket(tmp_dir, mdpocket_run1_args)

        # take the produced dx file and create the density grid in pdb format
        self._create_density_grid_file(tmp_dir, iso_value=self.iso_value)

        # We don't need all of this, but cluster pockets
        data, labels, _, _, _ = self._cluster_pockets(
            tmp_dir=tmp_dir,
            eps=self.eps,
            min_samples=self.min_samples,
            threshold=self.threshold,
            iso_value=self.iso_value,
        )

        # produces a load of pocket_n.pdb files based on the clusters identified by dbscan
        self._save_pocket_files(tmp_dir, data, labels)
        # run MD pocket the second time with a specified pocket - produce a pocket parameter fil
        # this should be done for each individual pocekt, in parallele
        # check whether the descriptors flag has been set
        # if _SFP.DESCRIPTORS in self.settings.additional.keys() and self.settings.additional[_SFP.DESCRIPTORS]:
        self._run_mdpocket_selected_pocket(tmp_dir)

        # save what's in the tmpdir, then remove tmpdir
        self._parse_output(tmp_dir)
        self._logger.log(
            f"Completed execution for {self.step_id} successfully", _LE.INFO
        )
        self._remove_temporary(tmp_dir)
