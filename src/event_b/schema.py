"""Machine-readable JSON Schemas for canonical Event-B entity records."""

from __future__ import annotations

import json
from pathlib import Path

from event_b.models import SCHEMAS


def entity_json_schema(entity: str) -> dict:
    schema = SCHEMAS[entity]
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": f"https://epicurus.local/schemas/{schema.version}/{entity}.json",
        "title": f"Epicurus Event-B {entity}",
        "type": "object",
        "additionalProperties": False,
        "required": list(schema.required),
        "properties": {
            column: ({"const": schema.version} if column == "schema_version" else {})
            for column in schema.columns
        },
    }


def write_schema_bundle(output_dir: str | Path) -> dict[str, str]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    written = {}
    for entity in sorted(SCHEMAS):
        path = output / f"{entity}.schema.json"
        path.write_text(json.dumps(entity_json_schema(entity), indent=2, sort_keys=True) + "\n")
        written[entity] = str(path)
    return written
