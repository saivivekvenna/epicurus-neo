"""Shared deterministic record helpers used by multiple real study adapters."""

from __future__ import annotations

from hashlib import sha256

import pandas as pd

from event_b.models import ReviewStatus, SCHEMAS, ValueOrigin


def stable_record_id(prefix: str, *parts: object) -> str:
    identity = "|".join(str(part) for part in parts)
    return f"{prefix}:" + sha256(identity.encode()).hexdigest()[:20]


def provenance_record(
    entity: str,
    entity_id: str,
    *,
    document: str,
    table: str,
    row: int | str,
    column: str,
    fragment: str,
    method: str = "deterministic_xlsx_adapter",
    origin: str = ValueOrigin.SOURCE_REPORTED.value,
) -> dict:
    provenance_id = stable_record_id("prov", entity, entity_id)
    return {
        "provenance_id": provenance_id,
        "entity_type": entity,
        "entity_id": entity_id,
        "field_name": "*",
        "source_document": document,
        "table": table,
        "row": str(row),
        "column": column,
        "source_fragment": fragment,
        "extraction_method": method,
        "extraction_confidence": 1.0,
        "value_origin": origin,
        "review_status": ReviewStatus.ACCEPTED.value,
    }


def entity_frame(entity: str, rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=SCHEMAS[entity].columns)
