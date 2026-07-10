"""Review issues are first-class records, never silently resolved."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path


@dataclass(frozen=True)
class ReviewIssue:
    issue_id: str
    entity_type: str
    entity_id: str
    code: str
    message: str
    severity: str = "ERROR"
    review_status: str = "NEEDS_REVIEW"
    candidate_record: dict | None = None
    conflicting_fields: tuple[str, ...] = ()
    source_evidence: tuple[dict, ...] = ()
    validation_rule: str | None = None
    suggested_resolution: str | None = None

    @classmethod
    def create(
        cls,
        entity_type: str,
        entity_id: str,
        code: str,
        message: str,
        *,
        candidate_record: dict | None = None,
        conflicting_fields: tuple[str, ...] = (),
        source_evidence: tuple[dict, ...] = (),
        suggested_resolution: str | None = None,
    ):
        identity = f"{entity_type}|{entity_id}|{code}|{message}"
        return cls(
            sha256(identity.encode()).hexdigest()[:20],
            entity_type,
            entity_id,
            code,
            message,
            candidate_record=candidate_record,
            conflicting_fields=conflicting_fields,
            source_evidence=source_evidence,
            validation_rule=code,
            suggested_resolution=suggested_resolution,
        )


def write_review_queue(issues: list[ReviewIssue], path: str | Path) -> None:
    rows = sorted((asdict(issue) for issue in issues), key=lambda row: row["issue_id"])
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in rows))


def read_review_queue(path: str | Path) -> list[ReviewIssue]:
    text = Path(path).read_text() if Path(path).exists() else ""
    rows = []
    for line in text.splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        row["conflicting_fields"] = tuple(row.get("conflicting_fields", ()))
        row["source_evidence"] = tuple(row.get("source_evidence", ()))
        rows.append(ReviewIssue(**row))
    return rows
