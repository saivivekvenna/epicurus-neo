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

    @classmethod
    def create(cls, entity_type: str, entity_id: str, code: str, message: str):
        identity = f"{entity_type}|{entity_id}|{code}|{message}"
        return cls(
            sha256(identity.encode()).hexdigest()[:20], entity_type, entity_id, code, message
        )


def write_review_queue(issues: list[ReviewIssue], path: str | Path) -> None:
    rows = sorted((asdict(issue) for issue in issues), key=lambda row: row["issue_id"])
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in rows))
