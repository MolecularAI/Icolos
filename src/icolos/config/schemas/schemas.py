import os
import json
from typing import Tuple

from icolos.utils.enums.general_utils_enums import JSONSchemasEnum
from icolos.utils.general.files_paths import attach_root_path

_JSE = JSONSchemasEnum()


def _load_schema(path: str) -> dict:
    try:
        with open(path, "r") as f:
            schema_data = f.read()
        schema = json.loads(schema_data)
    except:
        print(f"Could not load or parse schema file {path}, skipping.")
        schema = {}
    return schema


def _construct_absolute_path(sub_schema: str) -> str:
    rel_path = "src/icolos/config/schemas/"
    return attach_root_path(os.path.join(rel_path, "".join([sub_schema, ".json"])))


def construct_workflow_schema() -> Tuple[dict, str]:
    path = _construct_absolute_path(_JSE.WORKFLOW_SCHEMA)
    return _load_schema(path), os.path.dirname(path)
